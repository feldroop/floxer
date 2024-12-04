#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <seqan3/alphabet/cigar/cigar.hpp>

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

bool operator==(query_alignment const& lhs, query_alignment const& rhs);

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

    // returns the lowest number of errors among all alignments stored
    std::optional<size_t> best_num_errors() const;

    size_t size() const;

    // the other one is consumed (should be moved into this function)
    void merge_other_into_this(query_alignments other);
};

enum class alignment_mode {
    only_verify_existance, verify_and_return_alignment_with_cigar, verify_and_return_reduced_output
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

alignment_result align(
    std::span<const uint8_t> const reference,
    std::span<const uint8_t> const query,
    alignment_config const& config
);

} // namespace verification
