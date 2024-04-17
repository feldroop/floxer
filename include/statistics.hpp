#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <spdlog/spdlog.h>
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

    static inline const std::string query_lengths_name = "query lengths";
    static inline const std::string seed_lengths_name = "seed lengths";
    static inline const std::string anchors_per_seed_name = "anchors per seed";
    static inline const std::string anchors_per_query_name = "anchors per query";
    static inline const std::string alignments_per_query_name = "alignments per query";
    static inline const std::string alignments_edit_distance_name = "alignments edit distance";

    std::vector<histogram> histograms{
        histogram{large_values_log_scale, query_lengths_name},
        histogram{small_values_linear_scale, seed_lengths_name},
        histogram{large_values_log_scale, anchors_per_seed_name},
        histogram{large_values_log_scale, anchors_per_query_name},
        histogram{large_values_log_scale, alignments_per_query_name},
        histogram{small_values_log_scale, alignments_edit_distance_name}
    };

    histogram& histogram_by_name(std::string const& name);
    histogram const& histogram_by_name(std::string const& name) const;

    void insert_value_to(std::string const& target_name, size_t const value);

public:
    void add_query_length(size_t const value);

    void add_seed_length(size_t const value);

    void add_num_anchors_per_seed(size_t const value);

    void add_num_anchors_per_query(size_t const value);

    void add_num_alignments(size_t const value);

    void add_alignment_edit_distance(size_t const value);

    size_t num_queries() const;

    std::vector<std::string> format_histograms() const;

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
