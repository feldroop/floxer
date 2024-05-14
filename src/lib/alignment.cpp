#include <alignment.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <utility>

#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/adaptation/uint.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include <spdlog/spdlog.h>

#include <wavefront/wavefront_align.h>

namespace alignment {

query_alignments::query_alignments(size_t const num_references)
    : alignments_per_reference(num_references) {}

void query_alignments::insert(query_alignment const alignment, size_t const reference_id) {
    best_num_errors_ = std::min(
        best_num_errors_.value_or(std::numeric_limits<size_t>::max()),
        alignment.num_errors
    );

    alignments_per_reference[reference_id].emplace_back(std::move(alignment));
}

query_alignments::alignments_to_reference const& query_alignments::to_reference(
    size_t const reference_id
) const {
    return alignments_per_reference[reference_id];
}

query_alignments::alignments_to_reference& query_alignments::to_reference(size_t const reference_id) {
    return alignments_per_reference[reference_id];
}

std::optional<size_t> query_alignments::best_num_errors() const {
    return best_num_errors_;
}

size_t query_alignments::size() const {
    size_t size = 0;

    for (auto const& alignments_of_reference : alignments_per_reference) {
        size += alignments_of_reference.size();
    }

    return size;
}

// this must be global, because the aligner needs a default constructor for OpenMP private use
// and we wan't a new aligner to be constructed in every thread
static alignment_backend alignment_backend_global = alignment_backend::seqan3;

void set_alignment_backend_global(alignment_backend const backend) {
    alignment_backend_global = backend;
}

// this exists to allow seqan3 alignments with raw uint8_t's
template <seqan3::arithmetic score_type = int8_t>
class uint8_adaptation_scoring_scheme : public seqan3::scoring_scheme_base<
    uint8_adaptation_scoring_scheme<score_type>,
    uint8_t, score_type
> {
private:
    using base_t = seqan3::scoring_scheme_base<uint8_adaptation_scoring_scheme<score_type>, uint8_t, score_type>;
    friend base_t;

public:
    constexpr uint8_adaptation_scoring_scheme() noexcept
        : base_t(seqan3::match_score<score_type>{0}, seqan3::mismatch_score<score_type>{-1}) {};
};

aligner::aligner() : backend(alignment_backend_global) {
    switch (backend) {
        case alignment::alignment_backend::seqan3:
            // no extra setup yet for seqan3 aligner
            return;

        case alignment::alignment_backend::wfa2:
            setup_for_wfa2();
            return;

        default:
            throw std::runtime_error("(should be unreachable) internal bug in alignment backend choice");
    }
}

void aligner::setup_for_wfa2() {
    wavefront_aligner_attr_t attributes_only_score = wavefront_aligner_attr_default;

    attributes_only_score.distance_metric = distance_metric_t::edit;
    attributes_only_score.alignment_scope = alignment_scope_t::compute_score;

    // semi-global alignment, text_begin_free and text_end_free will be set for each alignment
    attributes_only_score.alignment_form.span = alignment_span_t::alignment_endsfree;
    attributes_only_score.alignment_form.pattern_begin_free = 0;
    attributes_only_score.alignment_form.pattern_end_free = 0;

    // we have to use this memory mode, because the longread alignment can become quite large
    // according to GitHub comments by the author, it's often just as fast as the other memory modes
    attributes_only_score.memory_mode = wavefront_memory_t::wavefront_memory_ultralow;

    // maybe add some heuristics later?
    attributes_only_score.heuristic.strategy = wf_heuristic_strategy::wf_heuristic_none;

    wavefront_aligner_attr_t attributes_full_alignment = attributes_only_score;
    attributes_only_score.alignment_scope = alignment_scope_t::compute_alignment;

    wf_aligner_only_score = wavefront_aligner_new(&attributes_only_score);
    wf_aligner_full_alignment = wavefront_aligner_new(&attributes_full_alignment);
}

alignment_result aligner::align(
    std::span<const uint8_t> const reference,
    std::span<const uint8_t> const query,
    alignment_config const& config
) {
    switch (backend) {
        case alignment::alignment_backend::seqan3:
            return align_seqan3(reference, query, config);

        case alignment::alignment_backend::wfa2:
            return align_wfa2(reference, query, config);

        default:
            throw std::runtime_error("(should be unreachable) internal bug in alignment backend choice");
    }
}

alignment_result aligner::align_seqan3(
    std::span<const uint8_t> const reference,
    std::span<const uint8_t> const query,
    alignment_config const& config
) {
    int32_t const min_score = -static_cast<int>(config.num_allowed_errors);
    auto aligner_config = seqan3::align_cfg::method_global{
        seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
        seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
        seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
        seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}
    }
    | seqan3::align_cfg::scoring_scheme{uint8_adaptation_scoring_scheme{}}; // gap scheme defaults to edit distance
    // | seqan3::align_cfg::min_score{min_score};

    if (config.mode == alignment_mode::only_verify_existance) {
        auto small_output_config = seqan3::align_cfg::output_score{};

        auto alignment_results = seqan3::align_pairwise(std::tie(reference, query), aligner_config | small_output_config);

        return alignment_result {
            .outcome = alignment_results.begin()->score() >= min_score ?
                alignment_outcome::alignment_exists :
                alignment_outcome::no_adequate_alignment_exists
        };
    }

    // when we need to do traceback, the matrix size might be a problem
    size_t const estimated_matrix_size = reference.size() * 2 * config.num_allowed_errors;
    size_t const very_big_matrix_size = 10'000'000'000;
    if (estimated_matrix_size > very_big_matrix_size) {
        spdlog::warn("Large alignment matrix of size {}", estimated_matrix_size);
    }

    auto full_output_config = seqan3::align_cfg::output_score{}
        | seqan3::align_cfg::output_begin_position{}
        | seqan3::align_cfg::output_alignment{};

    auto alignment_results = seqan3::align_pairwise(std::tie(reference, query), aligner_config | full_output_config);
    auto const alignment = *alignment_results.begin();

    if (alignment.score() < min_score) {
        return alignment_result { .outcome = alignment_outcome::no_adequate_alignment_exists };
    }

    assert(alignment.sequence2_begin_position() == 0);

    return alignment_result {
        .outcome = alignment_outcome::alignment_exists,
        .alignment = query_alignment {
            .start_in_reference = config.reference_span_offset + alignment.sequence1_begin_position(),
            .num_errors = static_cast<size_t>(std::abs(alignment.score())),
            .orientation = config.orientation,
            .cigar = seqan3::cigar_from_alignment(alignment.alignment(), {}, true)
        }
    };
}

alignment_result aligner::align_wfa2(
    std::span<const uint8_t> const reference,
    std::span<const uint8_t> const query,
    alignment_config const& config
) {
    // TODO

    // attributes_only_score.alignment_form.text_begin_free = 0;
    // attributes_only_score.alignment_form. = 0;

    return alignment_result { .outcome = alignment_outcome::no_adequate_alignment_exists };
}

aligner::~aligner() {
    wavefront_aligner_delete(wf_aligner_only_score);
    wavefront_aligner_delete(wf_aligner_full_alignment);
}


} // namespace alignment
