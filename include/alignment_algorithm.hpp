#pragma once

#include <alignment.hpp>
#include <concepts.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace alignment {

class VerifyingAligner {
public:
    // return whether an alignment between query and reference with at most the given number
    // of errors exists. if output_alignments is true, add all "useful" alignments to found_alignments.
    // useful means:
    //     - at most one for each position
    //     - must have at most num_allowed_errors errors
    //     - there is no locally better alignment
    //     - in case of multiple alignment per position apply the tiebreaking:
    //           mismatch > insertion to reference > deletion from from reference
    template<Sequence Seq>
    static bool align_query(
        Seq const& reference,
        Seq const& query,
        size_t const num_allowed_errors,
        bool const output_alignments,
        alignment_insertion_gatekeeper& alignment_gatekeeper
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

private:
    enum class trace_t {
        none, take_both_match, take_both_mismatch, only_query, only_reference 
    };

    using score_matrix_t = std::vector<std::vector<int>>;
    using traceback_matrix_t = std::vector<std::vector<trace_t>>;

    struct scoring_t {
        int const gap_score;
        int const mismatch_score;
        int const match_score;
    };

    static scoring_t edit_distance() {
        return scoring_t {
            .gap_score = -1,
            .mismatch_score = -1,
            .match_score = 0
        };
    }

    static std::tuple<score_matrix_t, traceback_matrix_t> initialize_matrices(
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

    template<Sequence Seq>
    static void fill_matrices(
        score_matrix_t& score_matrix,
        traceback_matrix_t& traceback_matrix,
        Seq const& reference,
        Seq const& query,
        scoring_t const& scoring
    ) {
        for (size_t i = 0; i < query.size(); ++i) {
            for (size_t j = 0; j < reference.size(); ++j) {
                int const take_only_reference_score = score_matrix[i + 1][j] + scoring.gap_score;      
                int const take_only_query_score = score_matrix[i][j + 1] + scoring.gap_score;      
                
                int score = score_matrix[i][j];
                trace_t trace;
                if (query[i] == reference[j]) {
                    score += scoring.match_score;
                    trace = trace_t::take_both_match;
                } else {
                    score += scoring.mismatch_score;
                    trace = trace_t::take_both_mismatch;
                }

                // the order of comparison guarantees the shortest alignment
                // regarding span over the reference:
                // first mismatch, then referernce insertion, then reference deletion
                if (take_only_reference_score > score) {
                    score = take_only_reference_score;
                    trace = trace_t::only_reference;
                }
                
                if (take_only_query_score > score) {
                    score = take_only_reference_score;
                    trace = trace_t::only_query;
                }

                score_matrix[i + 1][j + 1] = score;
                traceback_matrix[i + 1][j + 1] = trace;
            }
        }
    }

    static int find_best_score(score_matrix_t const& score_matrix) {
        auto const & last_row = score_matrix.back();
        int best_score = std::numeric_limits<int>::min();
        for (size_t i = 0; i < last_row.size(); ++i) {
            best_score = std::max(best_score, last_row[i]);
        }

        return best_score;
    }

    static query_alignment traceback(
        size_t const traceback_start_index,
        traceback_matrix_t const& traceback_matrix,
        size_t const num_errors
    ) {
        size_t i = traceback_matrix.size() - 1;
        size_t j = traceback_start_index;
        trace_t trace = traceback_matrix[i][j];
        cigar_sequence cigar;

        while (trace != trace_t::none) {
            switch (trace) {
                case trace_t::take_both_match:
                    cigar.add_operation(alignment_operation::match);
                    --i;
                    --j;
                    break;
                
                case trace_t::take_both_mismatch:
                    cigar.add_operation(alignment_operation::mismatch);
                    --i;
                    --j;
                    break;

                case trace_t::only_query:
                    cigar.add_operation(alignment_operation::deletion_from_reference);
                    --i;
                    break;

                case trace_t::only_reference:
                    cigar.add_operation(alignment_operation::insertion_to_reference);
                    --j;
                    break;
                
                default:
                    throw std::runtime_error("This should be unreachable.");
            }

            trace = traceback_matrix[i][j];
        }

        // no reversal of the cigar string here because the sequences are reversed

        return query_alignment {
            .start_in_reference = j, // will be transformed by insertion gatekeeper
            .end_in_reference = traceback_start_index, // will be transformed by insertion gatekeeper
            .reference_id = 0, // will be overridden by insertion gatekeeper
            .num_errors = num_errors,
            .score = -static_cast<int64_t>(num_errors),
            .is_reverse_complement = false, // will be overridden by insertion gatekeeper
            .cigar = cigar
        };
    }

    static void collect_and_write_alignments(
        score_matrix_t const& score_matrix,
        traceback_matrix_t const& traceback_matrix,
        size_t const num_allowed_errors,
        alignment_insertion_gatekeeper& alignment_gatekeeper
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

            alignment_gatekeeper.insert_alignment_if_its_useful(i, curr_num_errors, [&] () {
                return traceback(i, traceback_matrix, curr_num_errors);
            });
        }
    }
};

} // namespace alignment
