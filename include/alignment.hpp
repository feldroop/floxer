#pragma once

#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace alignment {

enum class alignment_operation {
    match, mismatch, deletion_from_reference, insertion_to_reference
};

char to_char(alignment_operation v);

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

    size_t length_in_reference() const;

    // this assumes that the alignments reference id is equivalent
    alignment_quality_comparison local_quality_comparison_versus(
        size_t const other_start_in_reference,
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
        size_t const reference_span_length,
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
    size_t const reference_span_length;
    bool const is_reverse_complement;
    fastq_query_alignments& useful_existing_alignments;
public:
    alignment_insertion_gatekeeper(
        size_t const reference_id_,
        size_t const reference_span_start_offset_,
        size_t const reference_span_length_,
        bool const is_reverse_complement,
        fastq_query_alignments& useful_existing_alignments_
    );

    // returns whether the alignment was added
    bool insert_alignment_if_its_useful(
        // this value should be the end position as seen by the alignment algorithm,
        // but it's actually the start position, because the sequence was reversed
        size_t const candidate_reference_span_reverse_end_position,
        size_t const candidate_num_erros,
        std::function<query_alignment()> const compute_alignment_to_reference_span
    );
};

} // namespace verification
