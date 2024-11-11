#include <alignment.hpp>

#include <algorithm>
#include <cassert>
#include <charconv>
#include <cmath>
#include <limits>
#include <system_error>
#include <utility>

#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/adaptation/uint.hpp>

#include <spdlog/fmt/fmt.h>
#include <spdlog/spdlog.h>

namespace alignment {

bool operator==(query_alignment const& lhs, query_alignment const& rhs) {
    if (lhs.cigar.size() != rhs.cigar.size()) {
        return false;
    }

    for (size_t i = 0; i < lhs.cigar.size(); ++i) {
        if (lhs.cigar.at(i) != rhs.cigar.at(i)) {
            return false;
        }
    }

    return lhs.num_errors == rhs.num_errors &&
        lhs.orientation == rhs.orientation &&
        lhs.start_in_reference == rhs.start_in_reference;
}

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

void query_alignments::merge_other_into_this(query_alignments other) {
    for (size_t reference_id = 0; reference_id < alignments_per_reference.size(); ++reference_id) {
        for (auto& alignment : other.to_reference(reference_id)) {
            insert(std::move(alignment), reference_id);
        }
    }
}

static constexpr uint64_t very_large_memory_usage = 10'000'000'000;

alignment_result align(
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
    | seqan3::align_cfg::edit_scheme
    | seqan3::align_cfg::min_score{min_score};

    if (config.mode == alignment_mode::only_verify_existance) {
        auto small_output_config = seqan3::align_cfg::output_score{};

        auto alignment_results = seqan3::align_pairwise(std::tie(reference, query), aligner_config | small_output_config);
        auto const alignment = *alignment_results.begin();

        // if no alignment could be found, the score is set to this value by the min_score configuration
        auto const infinite = std::numeric_limits<decltype(alignment.score())>::max();

        return alignment_result {
            .outcome = alignment.score() == infinite ?
                alignment_outcome::no_adequate_alignment_exists :
                alignment_outcome::alignment_exists
        };
    }

    size_t const reference_surplus_size = reference.size() >= query.size() ? reference.size() - query.size() : 0;
    size_t const estimated_band_size = 2 * config.num_allowed_errors + reference_surplus_size;
    size_t const estimated_matrix_size = reference.size() * estimated_band_size;
    if (estimated_matrix_size > very_large_memory_usage) {
        spdlog::warn("Large seqan3 alignment matrix of estimated size {}", estimated_matrix_size);
    }

    auto full_output_config = seqan3::align_cfg::output_score{}
        | seqan3::align_cfg::output_begin_position{}
        | seqan3::align_cfg::output_alignment{};

    auto alignment_results = seqan3::align_pairwise(std::tie(reference, query), aligner_config | full_output_config);
    auto const alignment = *alignment_results.begin();

    // if no alignment could be found, the score is set to this value by the min_score configuration
    auto const infinite = std::numeric_limits<decltype(alignment.score())>::max();

    if (alignment.score() == infinite) {
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

} // namespace alignment
