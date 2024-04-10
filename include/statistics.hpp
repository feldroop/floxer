#pragma once

#include <array>
#include <cstddef>
#include <string>

#include <spdlog/spdlog.h>
#include <spdlog/fmt/fmt.h>

namespace statistics {

class search_and_alignment_statistics {
    static constexpr size_t NUM_HISTOGRAM_BUCKETS = 9;
    static constexpr std::array<size_t, NUM_HISTOGRAM_BUCKETS> histogram_thresholds = {
        0, 1, 5, 10, 20, 100, 1'000, 10'000, 100'000
    };

    struct histogram {
        std::array<size_t, NUM_HISTOGRAM_BUCKETS> data{0};
        size_t num_values{0};

        void add_value(size_t const value) {
            ++num_values;

            for (size_t i = 0; i < NUM_HISTOGRAM_BUCKETS - 1; ++i) {
                if (value <= histogram_thresholds[i + 1]) {
                    ++data[i];
                    return;
                }
            }

            ++data[NUM_HISTOGRAM_BUCKETS - 1];
        }

        void merge_with(histogram const& other) {
            num_values += other.num_values;

            for (size_t i = 0; i < NUM_HISTOGRAM_BUCKETS; ++i) {
                data[i] += other.data[i];
            }
        }

        std::string format_to_string(std::string const& name) const {
            return fmt::format(
                "histogram for {} (total: {})\nthreshold:\t{}\noccurrences:\t{}",
                name,
                num_values,
                fmt::join(histogram_thresholds, "\t"),
                fmt::join(data, "\t")
            );
        }
    };

    histogram seed_hits_histogram{};
    histogram seed_hits_per_query_histogram{};
    histogram num_alignments_forward_histogram{};
    histogram num_alignments_revcomp_histogram{};

public:
    void add_num_seed_search_hits(size_t const value) {
        seed_hits_histogram.add_value(value);
    }

    void add_num_seed_search_hits_per_query(size_t const value) {
        seed_hits_per_query_histogram.add_value(value);
    }

    void add_num_alignments_forward(size_t const value) {
        num_alignments_forward_histogram.add_value(value);
    }

    void add_num_alignments_revcomp(size_t const value) {
        num_alignments_revcomp_histogram.add_value(value);
    }

    void print_all_histograms() const {
        spdlog::info(
            "{}", seed_hits_histogram.format_to_string("seed hits per seed")
        );
        spdlog::info(
            "{}", seed_hits_per_query_histogram.format_to_string("seed hits per query")
        );
        spdlog::info(
            "{}", num_alignments_forward_histogram.format_to_string("alignments forward")
        );
        spdlog::info(
            "{}", num_alignments_revcomp_histogram.format_to_string("alignments revcomp")
        );
    }

    friend search_and_alignment_statistics& combine_stats(
        search_and_alignment_statistics & inout,
        search_and_alignment_statistics const& other
    );
};

inline search_and_alignment_statistics& combine_stats(
    search_and_alignment_statistics & inout,
    search_and_alignment_statistics const& other
) {
    inout.seed_hits_histogram.merge_with(other.seed_hits_histogram);
    inout.seed_hits_per_query_histogram.merge_with(other.seed_hits_per_query_histogram);
    inout.num_alignments_forward_histogram.merge_with(other.num_alignments_forward_histogram);
    inout.num_alignments_revcomp_histogram.merge_with(other.num_alignments_revcomp_histogram);

    return inout;
}

} // namespace statistics
