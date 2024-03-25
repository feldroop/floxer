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
    match, mismatch, deletion_from_reference, insertion_to_reference
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
    size_t reference_id;
    size_t num_errors;
    int64_t score;
    bool is_reverse_complement;
    cigar_sequence cigar;

    size_t length_in_reference();

    // this assumes that the alignments reference id is equivalent
    alignment_quality_comparison local_quality_comparison_versus(
        size_t const other_end_in_reference,
        size_t const other_num_errors
    ) const;
};

// for testing
std::vector<alignment_operation> alignment_from_string(std::string const& s);

class alignment_insertion_gatekeeper;

// important invariant: there are only useful (locally optimal) alignments stored in this class
// only the alignment insertion gatekeeper should be used to insert alignments into this class
class fastq_query_alignments {
    using reference_alignments = std::map<size_t, query_alignment>;

    std::vector<reference_alignments> alignments_per_reference;

    std::optional<int64_t> primary_alignment_score = std::nullopt;
    std::optional<size_t> primary_alignment_reference_id = std::nullopt;
    std::optional<size_t> primary_alignment_end_position = std::nullopt;

    friend alignment_insertion_gatekeeper;

public:
    fastq_query_alignments(size_t const num_references);

    alignment_insertion_gatekeeper get_insertion_gatekeeper(
        size_t const reference_id,
        size_t const reference_span_start_offset,
        bool const is_reverse_complement
    );

    reference_alignments const& for_reference(size_t const reference_id) const;

    // primary alignment is the alignment with the highest score.
    // in case of ties, the one with the lowest position is chosen
    // and then the one with the lowest reference id
    bool is_primary_alignment(query_alignment const& alignment) const;

private:
    void update_primary_alignment(query_alignment const& alignment);
};

// makes sure that only useful alignments are inserted into the fastq_query_alignments class
// the other job of this class is to transform the indices of the computed alignment for a span
// the indices of the whole reference
class alignment_insertion_gatekeeper {
    size_t const reference_id;
    size_t const reference_span_start_offset;
    bool const is_reverse_complement;
    fastq_query_alignments& useful_existing_alignments;
public:
    alignment_insertion_gatekeeper(
        size_t const reference_id_,
        size_t const reference_span_start_offset_,
        bool const is_reverse_complement,
        fastq_query_alignments& useful_existing_alignments_
    );

    // returns whether the alignment was added
    bool insert_alignment_if_its_useful(
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
//     - there is no locally better alignment
//     - in case of multiple alignment per position apply the tiebreaking:
//           mismatch > insertion to reference > deletion from from reference
bool align_query(
    std::span<const uint8_t> reference,
    std::span<const uint8_t> query,
    size_t const num_allowed_errors,
    bool const output_alignments,
    alignment_insertion_gatekeeper& alignment_gatekeeper
);

} // namespace verification
