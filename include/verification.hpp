#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

namespace verification {

enum class alignment_variant {
    match, deletion, insertion
};

char format_as(alignment_variant v);

struct query_alignment {
    // half open range [start_in_reference, end_in_reference)
    size_t start_in_reference;
    size_t end_in_reference;
    size_t num_errors;

    // the terminolgy ov the alignment variants refers is from the perspective of the query
    // e.g. deletion -> deletion in the query
    std::vector<alignment_variant> alignment;

    size_t length_in_reference();
};

std::vector<alignment_variant> alignment_from_string(std::string const& s);

class full_reference_alignments {
    size_t const reference_span_start_offset;
    std::unordered_map<size_t, query_alignment>& alignments;
public:
    full_reference_alignments(
        size_t const reference_span_start_offset_,
        std::unordered_map<size_t, query_alignment>& alignments_
    ) : reference_span_start_offset{reference_span_start_offset_}, alignments{alignments_} 
    {}

    bool contains_equal_or_better_alignment_at_end_position(
        size_t const reference_span_end_position,
        size_t const new_alignment_num_erros
    ) const;

    void add(query_alignment && alignment_in_reference_span);
};

// return whether an alignment between query and reference with at most the given number
// of errors exists. if output_alignments is true, add all "useful" alignments to found_alignments.
// useful means:
//     - at most one for each position
//           (anchor is alignment end, TODO maybe make it the start)
//     - must have at most num_allowed_errors errors
//     - there is no directly neighboring alignment with strictly better score
//     - in case of multiple alignment per position apply the tiebreaking:
//           mismatch > insertion into query > deletion from query
bool align_query(
    std::span<const uint8_t> reference,
    std::span<const uint8_t> query,
    size_t const num_allowed_errors,
    bool const output_alignments,
    full_reference_alignments found_alignments
);

} // namespace verification
