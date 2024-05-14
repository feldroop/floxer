#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include <seqan3/alphabet/cigar/cigar.hpp>

// use the C API and not the C++ bindings of WFA2, because the C++ bindings
// don't allow retrieving the start offset of the alignment
extern "C" {
#include <wavefront/wavefront_align.h>
}

namespace alignment {

enum class query_orientation {
    forward, reverse_complement
};

struct query_alignment {
    size_t start_in_reference;
    size_t num_errors;
    query_orientation orientation;
    std::vector<seqan3::cigar> cigar;
};

// this class stores all of the alignments of one query to all references
class query_alignments {
    using alignments_to_reference = std::vector<query_alignment>;

    std::vector<alignments_to_reference> alignments_per_reference;

    std::optional<size_t> best_num_errors_ = std::nullopt;

public:
    query_alignments(size_t const num_references);

    void insert(query_alignment const alignment, size_t const reference_id);

    alignments_to_reference const& to_reference(size_t const reference_id) const;

    alignments_to_reference& to_reference(size_t const reference_id);

    std::optional<size_t> best_num_errors() const;

    size_t size() const;
};

enum class alignment_backend {
    seqan3, wfa2
};

void set_alignment_backend_global(alignment_backend const backend);

enum class alignment_mode {
    only_verify_existance, verify_and_return_alignment
};

struct alignment_config {
    size_t const reference_span_offset;
    size_t const num_allowed_errors;
    query_orientation const orientation;
    alignment_mode const mode;
};

enum class alignment_outcome {
    alignment_exists, no_adequate_alignment_exists
};

struct alignment_result {
    alignment_outcome outcome;
    std::optional<query_alignment> alignment = std::nullopt;
};

class aligner {
public:
    aligner();

    ~aligner();

    alignment_result align(
        std::span<const uint8_t> const reference,
        std::span<const uint8_t> const query,
        alignment_config const& config
    );

private:
    void setup_for_wfa2();

    alignment_result align_seqan3(
        std::span<const uint8_t> const reference,
        std::span<const uint8_t> const query,
        alignment_config const& config
    );

    alignment_result align_wfa2(
        std::span<const uint8_t> const reference,
        std::span<const uint8_t> const query,
        alignment_config const& config
    );

    void handle_wfa_status(const wavefront_align_status_t* const status);

    alignment_backend const backend;

    // only for wfa2 backend
    wavefront_aligner_t* wf_aligner_only_score = nullptr;
    wavefront_aligner_t* wf_aligner_full_alignment = nullptr;

    std::string wfa2_cigar_conversion_buffer{};
};

} // namespace verification
