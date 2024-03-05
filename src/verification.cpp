#include <verification.hpp>

#include <algorithm>
#include <cmath>
#include <compare>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <tuple>
#include <utility>

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

size_t query_alignment::length_in_reference() {
    return end_in_reference - start_in_reference;
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

constexpr scoring_t edit_distance() {
    return scoring_t {
        .gap_score = -1,
        .mismatch_score = -1,
        .match_score = 0
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

    auto & first_matrix_row = score_matrix.front();
    auto & first_traceback_row = traceback_matrix.front();
    for (size_t j = 0; j < first_matrix_row.size(); ++j) {
        // first_matrix_row[j] stays zero because of semiglobal alignment
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

std::vector<size_t> find_best_last_row_indices(
    score_matrix_t const& score_matrix,
    int const best_score
) {
    std::vector<size_t> best_last_row_indices{};
    auto const& last_row = score_matrix.back();
    for (size_t i = 0; i < last_row.size(); ++i) {
        if (last_row[i] == best_score) {
            best_last_row_indices.push_back(i);
        }
    }

    return best_last_row_indices;
}

query_alignment traceback(
    size_t const traceback_start_index,
    traceback_matrix_t const& traceback_matrix,
    size_t const best_num_errors
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
        .num_errors = best_num_errors,
        .alignment = alignment
    };
}

query_alignment tiebreak_alignment_choice(std::vector<query_alignment>& best_alignments) {
    query_alignment & shortest_alignment = best_alignments.front();
    for (size_t i = 1; i < best_alignments.size(); ++i) {
        size_t const current_shortest_length = shortest_alignment.length_in_reference();
        size_t const new_length = best_alignments[i].length_in_reference();
        auto const ordering = new_length <=> current_shortest_length;

        if (
            ordering == std::strong_ordering::less ||
            (
                ordering == std::strong_ordering::equal &&
                best_alignments[i].start_in_reference < shortest_alignment.start_in_reference
            )
        ) {
            shortest_alignment = best_alignments[i];
        }
    }

    return std::move(shortest_alignment);
}

std::optional<query_alignment> query_occurs(
    std::span<const uint8_t> const reference,
    std::span<const uint8_t> const query,
    size_t const num_allowed_errors
) {
    if (reference.empty() || query.empty()) {
        throw std::runtime_error("Empty sequences for verification alignment not allowed.");
    }

    scoring_t constexpr scoring = edit_distance();

    auto [score_matrix, traceback_matrix] = initialize_matrices(
        reference.size(), query.size(), scoring
    );

    fill_matrices(
        score_matrix, traceback_matrix, reference, query, scoring
    );

    int const best_score = find_best_score(score_matrix);

    size_t const best_num_errors = std::abs(best_score);
    if (best_num_errors > num_allowed_errors) {
        return std::nullopt;
    }

    std::vector<size_t> const best_last_row_indices = find_best_last_row_indices(score_matrix, best_score);

    std::vector<query_alignment> best_alignments{};
    for (size_t const last_row_index : best_last_row_indices) {
        auto alignment = traceback(
            last_row_index,
            traceback_matrix,
            best_num_errors
        );
        best_alignments.emplace_back(std::move(alignment));
    }

    return tiebreak_alignment_choice(best_alignments);
}

} // namespace verification
