#include <verification.hpp>

#include <algorithm>
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

bool full_reference_alignments::contains_equal_or_better_alignment_at_end_position(
    size_t const reference_span_end_position,
    size_t const new_alignment_num_erros
) const {
    size_t const end_in_full_reference =
        reference_span_start_offset + reference_span_end_position;
    
    auto const iter = alignments.find(end_in_full_reference);
    
    return iter != alignments.end() && iter->second.num_errors <= new_alignment_num_erros;
}

void full_reference_alignments::add(query_alignment && alignment_in_reference_span) {
    size_t const start_in_full_reference = 
        reference_span_start_offset + alignment_in_reference_span.start_in_reference;
    size_t const end_in_full_reference = 
        reference_span_start_offset + alignment_in_reference_span.end_in_reference;

    alignments[end_in_full_reference] = verification::query_alignment{
        .start_in_reference = start_in_full_reference,
        .end_in_reference = end_in_full_reference,
        .num_errors = alignment_in_reference_span.num_errors,
        .alignment = std::move(alignment_in_reference_span.alignment)
    };
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
    full_reference_alignments found_alignments
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

        if (
            found_alignments.contains_equal_or_better_alignment_at_end_position(i, curr_num_errors)
        ) {
            continue;
        }

        auto alignment = traceback(i, traceback_matrix, curr_num_errors);
        found_alignments.add(std::move(alignment));
    }
}

bool align_query(
    std::span<const uint8_t> const reference,
    std::span<const uint8_t> const query,
    size_t const num_allowed_errors,
    bool const output_alignments,
    full_reference_alignments found_alignments
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
            score_matrix, traceback_matrix, num_allowed_errors, found_alignments
        );
    }

    return true;
}

} // namespace verification
