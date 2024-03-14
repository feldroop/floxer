#pragma once

#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace verification {

enum class alignment_operation {
    match, mismatch, deletion, insertion
};

char format_as(alignment_operation v);

enum class alignment_quality_comparison {
    unrelated, equal, better, worse
};

class cigar_sequence {
    struct operation_block {
        alignment_operation operation;
        size_t count;
    };

    std::vector<operation_block> operation_blocks;

public:
    void add_operation(alignment_operation const operation);

    void reverse();

    std::string to_string() const;
};

struct query_alignment {
    // half open range [start_in_reference, end_in_reference)
    size_t start_in_reference;
    size_t end_in_reference;
    size_t num_errors;
    int64_t score;
    cigar_sequence cigar;

    size_t length_in_reference();

    alignment_quality_comparison local_quality_comparison_versus(
        size_t const other_end_in_reference,
        size_t const other_num_errors
    ) const;
};

std::string format_as(query_alignment const& alignment);

std::vector<alignment_operation> alignment_from_string(std::string const& s);

// important invariant: there are only useful (locally optimal) alignments stored in this class
// the other job of this class is to transform the indices of the computed alignment for a span
// the indices of the whole reference
class alignment_output_gatekeeper {
    size_t const reference_span_start_offset;
    std::map<size_t, query_alignment>& useful_existing_alignments;
public:
    alignment_output_gatekeeper(
        size_t const reference_span_start_offset_,
        std::map<size_t, query_alignment>& useful_existing_alignments_
    ) : reference_span_start_offset{reference_span_start_offset_},
        useful_existing_alignments{useful_existing_alignments_} 
    {}

    // returns whether the alignment was added
    bool add_alignment_if_its_useful(
        size_t const candidate_reference_span_end_position,
        size_t const candidate_num_erros,
        std::function<query_alignment()> const compute_alignment_to_reference_span
    );
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
    alignment_output_gatekeeper& alignment_gatekeeper
);

} // namespace verification
