#pragma once

#include <search.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <spdlog/fmt/fmt.h>

namespace statistics {

namespace internal {

std::vector<size_t> linear_range(size_t const num_steps, size_t const max);

} // namespace internal

class search_and_alignment_statistics {
public:
    struct histogram_config {
        std::vector<size_t> thresholds;
    };

    struct histogram_config_set {
        histogram_config small_values_linear_scale;
        histogram_config medium_values_linear_scale;
        histogram_config tiny_values_linear_scale;
        histogram_config practical_query_length_scale;
        histogram_config practical_anchor_scale;
        histogram_config kept_anchor_per_seed_scale;
        histogram_config edit_distance_scale;
        histogram_config practical_time_scale;
    };

private:
    struct count {
        std::string const name;
        size_t value = 0;

        std::string format_to_string_for_stdout() const;

        std::string format_as_toml() const;
    };

    struct histogram {
        histogram_config const config;
        std::string const name;

        std::vector<size_t> data{};
        size_t num_values{0};
        size_t min = std::numeric_limits<size_t>::max();
        double sum = 0.0;
        size_t max = std::numeric_limits<size_t>::min();

        histogram(histogram_config const& config_, std::string const& name_);

        void add_value(size_t const value);

        void merge_with(histogram const& other);

        std::string format_to_string_for_stdout() const;

        std::string format_as_toml() const;
    };

    static inline const std::string num_completely_excluded_queries_name = "completely excluded queries";

    std::vector<count> counts{
        count{num_completely_excluded_queries_name}
    };

    static inline const std::string query_lengths_name = "query lengths";
    static inline const std::string seed_lengths_name = "seed lengths";
    static inline const std::string errors_per_seed_name = "errors per seed";
    static inline const std::string seeds_per_query_name = "seeds per query";
    static inline const std::string anchors_per_seed_name = "anchors per non excluded seed";
    static inline const std::string kept_anchors_per_partly_excluded_seed_name = "kept anchors per partly excluded seed";
    static inline const std::string raw_anchors_per_excluded_seed_name = "raw anchors per fully excluded seed";
    static inline const std::string anchors_per_query_name = "anchors per query from non excluded seeds";
    static inline const std::string excluded_raw_anchors_per_query_name = "excluded raw anchors per query";
    static inline const std::string reference_span_sizes_aligned_inner_nodes_name = "reference span sizes aligned of inner nodes";
    static inline const std::string reference_span_sizes_aligned_root_name = "reference span sizes aligned of roots";
    static inline const std::string reference_span_sizes_avoided_root_name = "reference span sizes alignment avoided of roots";
    static inline const std::string alignments_per_query_name = "alignments per query";
    static inline const std::string alignments_edit_distance_name = "alignments edit distance";
    static inline const std::string milliseconds_spent_in_search_per_query_name = "milliseconds spent in search per query";
    static inline const std::string milliseconds_spent_in_verification_per_query_name = "milliseconds spent in verification per query";

    std::vector<histogram> histograms;

    count& count_by_name(std::string const& name);
    count const& count_by_name(std::string const& name) const;

    histogram& histogram_by_name(std::string const& name);
    histogram const& histogram_by_name(std::string const& name) const;

    void increment_count(std::string const& target_name);
    void insert_value_to(std::string const& target_name, size_t const value);

public:
    search_and_alignment_statistics(std::string_view const input_hint);

    void increment_num_completely_excluded_queries();

    void add_query_length(size_t const value);

    void add_seed_length(size_t const value);

    void add_num_errors_per_seed(size_t const value);

    void add_num_seeds_per_query(size_t const value);

    void add_statistics_for_seeds(
        std::vector<search::seed> const& forward_seeds,
        std::vector<search::seed> const& reverse_complement_seeds
    );

    void add_num_anchors_per_seed(size_t const value);

    void add_num_kept_anchors_per_partly_excluded_seed(size_t const value);

    void add_num_raw_anchors_per_excluded_seed(size_t const value);

    void add_num_anchors_per_query(size_t const value);

    void add_num_excluded_raw_anchors_per_query(size_t const value);

    void add_reference_span_size_aligned_inner_node(size_t const value);

    void add_reference_span_size_aligned_root(size_t const value);

    void add_reference_span_size_avoided_root(size_t const value);

    void add_num_alignments(size_t const value);

    void add_alignment_edit_distance(size_t const value);

    void add_milliseconds_spent_in_search_per_query(size_t const value);

    void add_milliseconds_spent_in_verification_per_query(size_t const value);

    void add_statistics_for_search_result(
        search::search_result const& forward_search_result,
        search::search_result const& reverse_complement_search_result
    );

    size_t num_queries() const;

    std::vector<std::string> format_statistics_for_stdout() const;

    std::string format_statistics_as_toml() const;

    // for merging local stats from threads into the global stats object
    void merge_other_into_this(search_and_alignment_statistics const& other);
};

} // namespace statistics
