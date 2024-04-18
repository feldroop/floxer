#pragma once

#include <search.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <spdlog/fmt/fmt.h>

namespace statistics {

class search_and_alignment_statistics {
    struct histogram_config {
        std::vector<size_t> thresholds;
    };

    static inline const histogram_config large_values_log_scale{
        .thresholds = { 0, 1, 5, 10, 20, 100, 1'000, 10'000, 100'000 }
    };

    static inline const histogram_config small_values_linear_scale{
        .thresholds = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70 }
    };

    static inline const histogram_config small_values_log_scale{
        .thresholds = { 0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000 }
    };

    struct count {
        std::string const name;
        size_t value = 0;

        std::string format_to_string() const;
    };

    struct histogram {
        histogram_config const config;
        std::string const name;

        std::vector<size_t> data{};
        size_t num_values{0};

        histogram(histogram_config const& config_, std::string const& name_);

        void add_value(size_t const value);

        void merge_with(histogram const& other);

        std::string format_to_string() const;
    };

    static inline const std::string num_completely_excluded_queries_name = "completely excluded queries";

    std::vector<count> counts{
        count{num_completely_excluded_queries_name}
    };

    static inline const std::string query_lengths_name = "query lengths";
    static inline const std::string seed_lengths_name = "seed lengths";
    static inline const std::string anchors_per_seed_name = "anchors per (non-excluded) seed";
    static inline const std::string raw_anchors_per_excluded_seed_name = "(raw) anchors per excluded seed";
    static inline const std::string anchors_per_query_name = "anchors per query (from non-excluded seeds)";
    static inline const std::string excluded_raw_anchors_per_query_name = "excluded (raw) anchors per query";
    static inline const std::string alignments_per_query_name = "alignments per query";
    static inline const std::string alignments_edit_distance_name = "alignments edit distance";

    std::vector<histogram> histograms{
        histogram{large_values_log_scale, query_lengths_name},
        histogram{small_values_linear_scale, seed_lengths_name},
        histogram{large_values_log_scale, anchors_per_seed_name},
        histogram{large_values_log_scale, raw_anchors_per_excluded_seed_name},
        histogram{large_values_log_scale, anchors_per_query_name},
        histogram{large_values_log_scale, excluded_raw_anchors_per_query_name},
        histogram{large_values_log_scale, alignments_per_query_name},
        histogram{small_values_log_scale, alignments_edit_distance_name}
    };

    count& count_by_name(std::string const& name);
    count const& count_by_name(std::string const& name) const;

    histogram& histogram_by_name(std::string const& name);
    histogram const& histogram_by_name(std::string const& name) const;

    void increment_count(std::string const& target_name);
    void insert_value_to(std::string const& target_name, size_t const value);

public:
    void increment_num_completely_excluded_queries();

    void add_query_length(size_t const value);

    void add_seed_length(size_t const value);

    void add_seed_lengths(std::vector<search::seed> const& seeds);

    void add_num_anchors_per_seed(size_t const value);

    void add_num_raw_anchors_per_excluded_seed(size_t const value);

    void add_num_anchors_per_query(size_t const value);
    
    void add_num_excluded_raw_anchors_per_query(size_t const value);

    void add_num_alignments(size_t const value);

    void add_alignment_edit_distance(size_t const value);

    void add_statistics_for_search_result(search::search_result const& search_result);

    size_t num_queries() const;

    std::vector<std::string> format_statistics() const;

    friend search_and_alignment_statistics& combine_stats(
        search_and_alignment_statistics & inout,
        search_and_alignment_statistics const& other
    );
};

search_and_alignment_statistics& combine_stats(
    search_and_alignment_statistics & inout,
    search_and_alignment_statistics const& other
);

} // namespace statistics
