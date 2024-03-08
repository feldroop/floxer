#include <verification.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <tuple>
#include <utility>

#include <fmt/core.h>
#include <fmt/ranges.h>

namespace verification {

enum class trace_t {
    none, take_both, only_query, only_reference 
};

struct scoring_t {
    int const gap_score;
    int const mismatch_score;
    int const match_score;
};

using score_matrix_t = std::vector<std::vector<int>>;
using traceback_matrix_t = std::vector<std::vector<trace_t>>;

char format_as(alignment_variant v) {
    switch (v) {
        case alignment_variant::match:
            return 'M';
        
        case alignment_variant::deletion:
            return 'D';

        case alignment_variant::insertion:
            return 'I';
        
        default:
            return 'X';
    }
}

size_t query_alignment::length_in_reference() {
    return end_in_reference - start_in_reference;
}

alignment_quality_comparison query_alignment::local_quality_comparison_versus(
    size_t const other_end_in_reference,
    size_t const other_num_errors
) const {
    size_t const distance_in_reference = end_in_reference >= other_end_in_reference ? 
        end_in_reference - other_end_in_reference :
        other_end_in_reference - end_in_reference;
    
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

std::string format_as(query_alignment const& alignment) {
    return fmt::format(
        "[{}, {}), {} errors, {}",
        alignment.start_in_reference,
        alignment.end_in_reference,
        alignment.num_errors,
        alignment.alignment
    );
}

std::vector<alignment_variant> alignment_from_string(std::string const& s) {
    std::vector<alignment_variant> alignment{};

    for (char const c : s) {
        switch (c) {
            case 'M':
                alignment.push_back(alignment_variant::match);
                break;
            
            case 'I':
                alignment.push_back(alignment_variant::insertion);
                break;

            case 'D':
                alignment.push_back(alignment_variant::deletion);
                break;

            default:
                throw std::runtime_error("Unexpected character in string");
        }
    }

    return alignment;
}

scoring_t edit_distance() {
    return scoring_t {
        .gap_score = -1,
        .mismatch_score = -1,
        .match_score = 0
    };
}

bool alignment_output_gatekeeper::add_alignment_if_its_useful(
    size_t const candidate_reference_span_end_position,
    size_t const candidate_num_erros,
    std::function<query_alignment()> const compute_alignment_to_reference_span
) {
    size_t const candidate_end_in_full_reference =
        reference_span_start_offset + candidate_reference_span_end_position;
    
    std::optional<size_t> worse_to_the_right_key = std::nullopt;

    // first check to the right
    auto iter = useful_existing_alignments.lower_bound(candidate_end_in_full_reference);
    if (iter != useful_existing_alignments.end()) {
        // this alignment might be at the exact same position or to the right of the candidate
        auto const& existing_alignment = iter->second;
        auto const existing_alignment_quality_ordering = existing_alignment.local_quality_comparison_versus(
            candidate_end_in_full_reference,
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
    if (iter != useful_existing_alignments.begin()) {
        --iter;
        // this alignment is to the left of the candidate
        auto const& existing_alignment = iter->second;
        auto const existing_alignment_quality_ordering = existing_alignment.local_quality_comparison_versus(
            candidate_end_in_full_reference,
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
        useful_existing_alignments.erase(iter);
    }

    if (worse_to_the_right_key.has_value()) {
        useful_existing_alignments.erase(worse_to_the_right_key.value());
    }

    auto reference_span_alignment = compute_alignment_to_reference_span();
    useful_existing_alignments.emplace(
        candidate_end_in_full_reference,
        query_alignment {
            .start_in_reference = reference_span_start_offset + reference_span_alignment.start_in_reference,
            .end_in_reference = reference_span_start_offset + reference_span_alignment.end_in_reference,
            .num_errors = reference_span_alignment.num_errors,
            .alignment = std::move(reference_span_alignment.alignment)
        }
    );

    return true;
}

std::tuple<score_matrix_t, traceback_matrix_t> initialize_matrices(
    size_t const reference_size,
    size_t const query_size,
    scoring_t const& scoring
) {
    score_matrix_t score_matrix(
        query_size + 1, 
        std::vector<int>(reference_size + 1, 0)
    );

    traceback_matrix_t traceback_matrix(
        query_size + 1, 
        std::vector<trace_t>(reference_size + 1, trace_t::none)
    );

    for (size_t i = 0; i < score_matrix.size(); ++i) {
        score_matrix[i][0] = i * scoring.gap_score;
        traceback_matrix[i][0] = trace_t::only_query;
    }

    auto & first_traceback_row = traceback_matrix.front();
    for (size_t j = 0; j < first_traceback_row.size(); ++j) {
        // first_score_matrix_row[j] stays zero because of semiglobal alignment
        first_traceback_row[0] = trace_t::none;
    }

    return std::make_tuple(std::move(score_matrix), std::move(traceback_matrix));
};

void fill_matrices(
    score_matrix_t& score_matrix,
    traceback_matrix_t& traceback_matrix,
    std::span<const uint8_t> const reference,
    std::span<const uint8_t> const query,
    scoring_t const& scoring
) {
    for (size_t i = 0; i < query.size(); ++i) {
        for (size_t j = 0; j < reference.size(); ++j) {
            int const take_only_reference_score = score_matrix[i + 1][j] + scoring.gap_score;      
            int const take_only_query_score = score_matrix[i][j + 1] + scoring.gap_score;      
            int const take_both_score = score_matrix[i][j] +
                (query[i] == reference[j] ? scoring.match_score : scoring.mismatch_score);

            // the order of comparison guarantees the shortest alignment
            // regarding span over the reference:
            // first insertion in query, then match/mismatch, then delete
            int score = take_only_query_score;
            trace_t trace = trace_t::only_query;

            if (take_both_score > score) {
                score = take_both_score;
                trace = trace_t::take_both;
            }
            
            if (take_only_reference_score > score) {
                score = take_only_reference_score;
                trace = trace_t::only_reference;
            }

            score_matrix[i + 1][j + 1] = score;
            traceback_matrix[i + 1][j + 1] = trace;
        }
    }
};

int find_best_score(score_matrix_t const& score_matrix) {
    auto const & last_row = score_matrix.back();
    int best_score = std::numeric_limits<int>::min();
    for (size_t i = 0; i < last_row.size(); ++i) {
        best_score = std::max(best_score, last_row[i]);
    }

    return best_score;
}

query_alignment traceback(
    size_t const traceback_start_index,
    traceback_matrix_t const& traceback_matrix,
    size_t const num_errors
) {
    size_t i = traceback_matrix.size() - 1;
    size_t j = traceback_start_index;
    trace_t trace = traceback_matrix[i][j];
    std::vector<alignment_variant> alignment{};

    while (trace != trace_t::none) {
        switch (trace) {
            case trace_t::take_both:
                alignment.push_back(alignment_variant::match);
                --i;
                --j;
                break;

            case trace_t::only_query:
                alignment.push_back(alignment_variant::insertion);
                --i;
                break;

            case trace_t::only_reference:
                alignment.push_back(alignment_variant::deletion);
                --j;
                break;
            
            default:
                throw std::runtime_error("This should be unreachable.");
        }

        trace = traceback_matrix[i][j];
    }

    std::ranges::reverse(alignment);

    return query_alignment {
        .start_in_reference = j,
        .end_in_reference = traceback_start_index,
        .num_errors = num_errors,
        .alignment = alignment
    };
}

void collect_and_write_alignments(
    score_matrix_t const& score_matrix,
    traceback_matrix_t const& traceback_matrix,
    size_t const num_allowed_errors,
    alignment_output_gatekeeper& alignment_gatekeeper
) {
    auto const& last_row_scores = score_matrix.back();

    for (size_t i = 0; i < last_row_scores.size(); ++i) {
        size_t const curr_num_errors = std::abs(last_row_scores[i]);
        if (curr_num_errors > num_allowed_errors) {
            continue;
        }

        size_t const left_neighbor_index = i == 0 ? i : i - 1;
        size_t const right_neighbor_index = i == last_row_scores.size() - 1 ? i : i + 1;

        if (
            last_row_scores[i] < last_row_scores[left_neighbor_index] ||
            last_row_scores[i] < last_row_scores[right_neighbor_index]
        ) {
            continue;
        }

        alignment_gatekeeper.add_alignment_if_its_useful(i, curr_num_errors, [&] () {
            return traceback(i, traceback_matrix, curr_num_errors);
        });
    }
}

bool align_query(
    std::span<const uint8_t> const reference,
    std::span<const uint8_t> const query,
    size_t const num_allowed_errors,
    bool const output_alignments,
    alignment_output_gatekeeper& alignment_gatekeeper
) {
    if (reference.empty() || query.empty()) {
        throw std::runtime_error("Empty sequences for verification alignment not allowed.");
    }

    scoring_t const scoring = edit_distance();

    auto [score_matrix, traceback_matrix] = initialize_matrices(
        reference.size(), query.size(), scoring
    );

    fill_matrices(
        score_matrix, traceback_matrix, reference, query, scoring
    );

    int const best_score = find_best_score(score_matrix);

    size_t const best_num_errors = std::abs(best_score);
    if (best_num_errors > num_allowed_errors) {
        return false;
    }

    if (output_alignments) {
        collect_and_write_alignments(
            score_matrix, traceback_matrix, num_allowed_errors, alignment_gatekeeper
        );
    }

    return true;
}

} // namespace verification
