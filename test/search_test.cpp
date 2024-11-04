#include <fmindex.hpp>
#include <search.hpp>

#include <gtest/gtest.h>

TEST(search, search_seeds) {
    search::search_config config {
        .max_num_anchors = 10,
        .anchor_group_order = search::anchor_group_order_t::hybrid
    };

    search::internal::search_scheme_cache scheme_cache;

    std::vector<std::vector<uint8_t>> const references {
        { 1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4 },
        { 1,2,3,4,1,2,3,4 }
    };
    size_t const num_reference_sequences = references.size();

    size_t const suffix_array_sampling_rate  = 4;
    size_t const num_threads = 4;
    fmindex index(
        references,
        suffix_array_sampling_rate,
        num_threads
    );

    search::searcher searcher {
        .index = index,
        .num_reference_sequences = num_reference_sequences,
        .config = config
    };

    std::vector<uint8_t> const query {
        1,1,1,1,1,1, // matches exactly
        2,2,2,3,2,2, // matches uniquely with 1 mismatch
        1,2,3,1,2,3, // matches uniquely with 1 deletion
        4,3,2,1,4,2  // does not match
    };
    std::span<const uint8_t> query_span(query);

    std::vector<search::seed> seeds{
        search::seed {
            .sequence = query_span.subspan(0,6),
            .num_errors = 0,
            .query_position = 0
        },
        search::seed {
            .sequence = query_span.subspan(6,6),
            .num_errors = 1,
            .query_position = 6
        },
        search::seed {
            .sequence = query_span.subspan(12,6),
            .num_errors = 1,
            .query_position = 12
        },
        search::seed {
            .sequence = query_span.subspan(18,6),
            .num_errors = 0,
            .query_position = 18
        }
    };

    auto const result = searcher.search_seeds(seeds);

    EXPECT_EQ(result.num_fully_excluded_seeds, 0);

    using anchors_of_seed = search::search_result::anchors_of_seed;
    std::vector<anchors_of_seed> const expected_anchors_by_seed {
        anchors_of_seed {
            .status = search::seed_status::not_excluded,
            .num_kept_useful_anchors = 1,
            .num_excluded_raw_anchors = 0,
            .anchors_by_reference = std::vector<search::anchors> {
                search::anchors {
                    search::anchor_t {
                        .reference_position = 0,
                        .num_errors = 0
                    }
                },
                {}
            }
        },
        anchors_of_seed {
            .status = search::seed_status::not_excluded,
            .num_kept_useful_anchors = 1,
            .num_excluded_raw_anchors = 0,
            .anchors_by_reference = std::vector<search::anchors> {
                search::anchors {
                    search::anchor_t {
                        .reference_position = 6,
                        .num_errors = 1
                    }
                },
                {}
            }
        },
        anchors_of_seed {
            .status = search::seed_status::not_excluded,
            .num_kept_useful_anchors = 1,
            .num_excluded_raw_anchors = 0,
            .anchors_by_reference = std::vector<search::anchors> {
                {},
                search::anchors {
                    search::anchor_t {
                        .reference_position = 0,
                        .num_errors = 1
                    }
                }
            }
        },
        anchors_of_seed {
            .status = search::seed_status::not_excluded,
            .num_kept_useful_anchors = 0,
            .num_excluded_raw_anchors = 0,
            .anchors_by_reference = std::vector<search::anchors> {
                {},
                {}
            }
        }
    };
}

TEST(search, erase_useless_anchors) {
    search::anchor_t const useful_anchor1 {
        .reference_position = 100,
        .num_errors = 0
    };

    search::anchor_t const useful_anchor2 {
        .reference_position = 120,
        .num_errors = 0
    };

    std::vector<search::anchors> anchors {{
        search::anchor_t {
            .reference_position = 95,
            .num_errors = 5
        },
        search::anchor_t {
            .reference_position = 97,
            .num_errors = 3
        },
        useful_anchor1,
        search::anchor_t {
            .reference_position = 110,
            .num_errors = 10
        },
        useful_anchor2
    }};


    search::internal::erase_useless_anchors(anchors);

    std::vector<search::anchors> const expected_anchors{{
        useful_anchor1, useful_anchor2
    }};

    EXPECT_EQ(anchors, expected_anchors);
}
