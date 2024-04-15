#include <alignment.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>

#include <spdlog/spdlog.h>

namespace alignment {

char to_char(alignment_operation v) {
    switch (v) {
        case alignment_operation::match:
            return '=';
        
        case alignment_operation::mismatch:
            return 'X';

        case alignment_operation::deletion_from_reference:
            return 'D';

        case alignment_operation::insertion_to_reference:
            return 'I';
        
        default:
            throw std::runtime_error("Unexpected alignment operation");
    }
}

void cigar_sequence::add_operation(alignment_operation const operation) {
    if (operation_blocks.empty() || operation_blocks.back().operation != operation) {
        operation_blocks.emplace_back(operation, 1);
    } else {
        operation_blocks.back().count += 1;
    }
}

void cigar_sequence::reverse() {
    std::ranges::reverse(operation_blocks);
}

size_t cigar_sequence::num_operation_blocks() const {
    return operation_blocks.size();
}

std::string cigar_sequence::to_string() const {
    std::stringstream stream{};

    for (auto const& block : operation_blocks) {
        stream << block.count << to_char(block.operation);
    }

    return stream.str();
}

size_t query_alignment::length_in_reference() const {
    return end_in_reference - start_in_reference;
}

alignment_quality_comparison query_alignment::local_quality_comparison_versus(
    size_t const other_start_in_reference,
    size_t const other_num_errors
) const {
    size_t const distance_in_reference = start_in_reference >= other_start_in_reference ? 
        start_in_reference - other_start_in_reference :
        other_start_in_reference - start_in_reference;
    
    size_t num_errors_difference;
    alignment_quality_comparison potential_comparison;

    if (num_errors > other_num_errors) {
        num_errors_difference = num_errors - other_num_errors;
        potential_comparison = alignment_quality_comparison::worse;
    } else if (num_errors == other_num_errors) {
        num_errors_difference = 0;
        potential_comparison = alignment_quality_comparison::equal;
    } else {
        num_errors_difference = other_num_errors - num_errors;
        potential_comparison = alignment_quality_comparison::better;
    }
    
    if (distance_in_reference > num_errors_difference) {
        return alignment_quality_comparison::unrelated;
    }

    return potential_comparison;
}

query_alignments::query_alignments(size_t const num_references)
    : alignments_per_reference(num_references) {}

alignment_insertion_gatekeeper query_alignments::get_insertion_gatekeeper(
    size_t const reference_id,
    size_t const reference_span_start_offset,
    size_t const reference_span_length,
    bool const is_reverse_complement
) {
    return alignment_insertion_gatekeeper(
        reference_id,
        reference_span_start_offset,
        reference_span_length,
        is_reverse_complement,
        *this
    );
}

bool query_alignments::is_primary_alignment(query_alignment const& alignment) const {
    if (!primary_alignment_end_position.has_value() || !primary_alignment_reference_id.has_value()) {
        return false;
    }

    return alignment.reference_id == primary_alignment_reference_id.value() &&
        alignment.end_in_reference == primary_alignment_end_position.value();
}

void query_alignments::update_primary_alignment(query_alignment const& new_alignment) {
    if (
        !primary_alignment_score.has_value() ||
        !primary_alignment_end_position.has_value() ||
        !primary_alignment_reference_id.has_value() ||
        primary_alignment_score.value() < new_alignment.score ||
        (
            primary_alignment_score.value() == new_alignment.score &&
            primary_alignment_end_position.value() > new_alignment.end_in_reference
        ) ||
        (
            primary_alignment_score.value() == new_alignment.score &&
            primary_alignment_end_position.value() == new_alignment.end_in_reference &&
            primary_alignment_reference_id.value() > new_alignment.reference_id
        )
    ) {
        primary_alignment_score = new_alignment.score;
        primary_alignment_end_position = new_alignment.end_in_reference;
        primary_alignment_reference_id = new_alignment.reference_id;
    }
}

size_t query_alignments::size() const {
    size_t size = 0;

    for (auto const& alignments_of_reference : alignments_per_reference) {
        size += alignments_of_reference.size();
    }

    return size;
}

query_alignments::alignments_to_reference const& query_alignments::to_reference(
    size_t const reference_id
) const {
    return alignments_per_reference[reference_id];
}

alignment_insertion_gatekeeper::alignment_insertion_gatekeeper(
    size_t const reference_id_,
    size_t const reference_span_start_offset_,
    size_t const reference_span_length_,
    bool const is_reverse_complement_,
    query_alignments& useful_existing_alignments_
) : reference_id{reference_id_},
    reference_span_start_offset{reference_span_start_offset_},
    reference_span_length{reference_span_length_},
    is_reverse_complement{is_reverse_complement_},
    useful_existing_alignments{useful_existing_alignments_} {}

bool alignment_insertion_gatekeeper::insert_alignment_if_its_useful(
    size_t const candidate_reference_span_reverse_end_position,
    size_t const candidate_num_erros,
    std::function<query_alignment()> const compute_alignment_to_reference_span
) {
    auto& this_reference_existing_alignments = 
        useful_existing_alignments.alignments_per_reference[this->reference_id];

    // first undo sequence reversal
    // in this term a - 1 + 1 cancelled out. The - 1 is from the backtransformation
    // and the + 1 is from the convention that end positions point to the first index
    // AFTER the sequence and start positions point to the first index INSIDE the sequence
    size_t const candidate_reference_span_start_position = 
        reference_span_length - candidate_reference_span_reverse_end_position;

    // then backtransform to whole sequence
    size_t const candidate_start_in_full_reference =
        reference_span_start_offset + candidate_reference_span_start_position;
    
    std::optional<size_t> worse_to_the_right_key = std::nullopt;

    // first check to the right
    auto iter = this_reference_existing_alignments.lower_bound(candidate_start_in_full_reference);
    if (iter != this_reference_existing_alignments.end()) {
        // this alignment might be at the exact same position or to the right of the candidate
        auto const& existing_alignment = iter->second;
        auto const existing_alignment_quality_ordering = existing_alignment.local_quality_comparison_versus(
            candidate_start_in_full_reference,
            candidate_num_erros
        );
        
        switch (existing_alignment_quality_ordering) {
            case alignment_quality_comparison::unrelated:
                // new alignment could be added later
                break;
            case alignment_quality_comparison::equal:
                return false;
            case alignment_quality_comparison::better:
                return false;
            case alignment_quality_comparison::worse:
                // new alignment will be added later and existing one will be deleted
                worse_to_the_right_key = iter->first;
                break;
            default:
                throw std::runtime_error("alignment quality comparison, this should be unreachable");
        }
    }

    bool found_worse_to_the_left = false;
    // now check the left
    if (iter != this_reference_existing_alignments.begin()) {
        --iter;
        // this alignment is to the left of the candidate
        auto const& existing_alignment = iter->second;
        auto const existing_alignment_quality_ordering = existing_alignment.local_quality_comparison_versus(
            candidate_start_in_full_reference,
            candidate_num_erros
        );

        switch (existing_alignment_quality_ordering) {
            case alignment_quality_comparison::unrelated:
                // new alignment will be added later
                break;
            case alignment_quality_comparison::better:
                // if this assert triggers, the invariant of the class was violated
                assert(!worse_to_the_right_key.has_value());
                return false;
            case alignment_quality_comparison::worse:
                // new alignment will be added later and existing one will be deleted
                found_worse_to_the_left = true;
                break;
            // case alignment_quality_comparison::equal: <- shouln't be possible
            default:
                throw std::runtime_error("alignment quality comparison, this should be unreachable");
        }
    }

    if (found_worse_to_the_left) {
        this_reference_existing_alignments.erase(iter);
    }

    if (worse_to_the_right_key.has_value()) {
        this_reference_existing_alignments.erase(worse_to_the_right_key.value());
    }

    auto reference_span_alignment = compute_alignment_to_reference_span();

    auto [inserted_iter, _] = this_reference_existing_alignments.emplace(
        candidate_start_in_full_reference,
        query_alignment {
            .start_in_reference = candidate_start_in_full_reference,
            .end_in_reference = candidate_start_in_full_reference + reference_span_alignment.length_in_reference(),
            .reference_id = this->reference_id,
            .num_errors = reference_span_alignment.num_errors,
            .score = reference_span_alignment.score,
            .is_reverse_complement = this->is_reverse_complement,
            .cigar = std::move(reference_span_alignment.cigar)
        }
    );

    useful_existing_alignments.update_primary_alignment(inserted_iter->second);

    return true;
}

} // namespace alignment
