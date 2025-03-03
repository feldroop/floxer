#include <pex.hpp>
#include <search.hpp>

#include <gtest/gtest.h>
#include <iostream>

TEST(pex, generate_seeds_from_recursive) {
    // simple AAACCCGGGTTT example

    // first, original pex with 0 leaf errors
    pex::pex_tree_config const config(
        12,
        3,
        0,
        pex::pex_tree_build_strategy::recursive
    );

    auto const tree = pex::pex_tree(config);

    std::vector<uint8_t> const query{ 0,0,0,1,1,1,2,2,2,3,3,3 };
    std::span<const uint8_t> query_span(query);
    auto const seeds = tree.generate_seeds(query_span, 1);

    std::vector<search::seed> expected_seeds{
        search::seed {
            .sequence = query_span.subspan(0,3),
            .num_errors = 0,
            .query_position = 0,
            .pex_leaf_index = 0
        },
        search::seed {
            .sequence = query_span.subspan(3,3),
            .num_errors = 0,
            .query_position = 3,
            .pex_leaf_index = 1
        },
        search::seed {
            .sequence = query_span.subspan(6,3),
            .num_errors = 0,
            .query_position = 6,
            .pex_leaf_index = 2
        },
        search::seed {
            .sequence = query_span.subspan(9,3),
            .num_errors = 0,
            .query_position = 9,
            .pex_leaf_index = 3
        },
    };
    EXPECT_EQ(seeds, expected_seeds);

    // now adjusted version with 1 leaf error
    pex::pex_tree_config const adjusted1_config(
        12,
        3,
        1,
        pex::pex_tree_build_strategy::recursive
    );
    auto const adjusted1_tree = pex::pex_tree(adjusted1_config);
    auto const adjusted1_seeds = adjusted1_tree.generate_seeds(query_span, 1);

    std::vector<search::seed> expected_seeds_adjusted {
        search::seed {
            .sequence = query_span.subspan(0,6),
            .num_errors = 1,
            .query_position = 0,
            .pex_leaf_index = 0
        },
        search::seed {
            .sequence = query_span.subspan(6,6),
            .num_errors = 1,
            .query_position = 6,
            .pex_leaf_index = 1
        },
    };
    EXPECT_EQ(adjusted1_seeds, expected_seeds_adjusted);

    // now adjusted version with 2 leaf errors, nothing should change
    pex::pex_tree_config const adjusted2_config(
        12,
        3,
        2,
        pex::pex_tree_build_strategy::recursive
    );
    auto const adjusted2_tree = pex::pex_tree(adjusted2_config);
    auto const adjusted2_seeds = adjusted2_tree.generate_seeds(query_span, 1);

    EXPECT_EQ(adjusted2_seeds, expected_seeds_adjusted);

}

TEST(pex, generate_seeds_from_bottom_up) {
    pex::pex_tree_config const config(
        30,
        14,
        2,
        pex::pex_tree_build_strategy::bottom_up
    );

    auto const tree = pex::pex_tree(config);

    std::vector<uint8_t> const query{
        0,0,0,1,1,1,2,2,2,3,
        3,3,0,0,0,1,1,1,2,2,
        2,3,3,3,0,0,0,1,1,1
    };
    std::span<const uint8_t> query_span(query);
    auto const seeds = tree.generate_seeds(query_span, 1);

    std::vector<search::seed> expected_seeds{
        search::seed {
            .sequence = query_span.subspan(0,6),
            .num_errors = 2,
            .query_position = 0,
            .pex_leaf_index = 0
        },
        search::seed {
            .sequence = query_span.subspan(6,6),
            .num_errors = 2,
            .query_position = 6,
            .pex_leaf_index = 1
        },
        search::seed {
            .sequence = query_span.subspan(12,6),
            .num_errors = 2,
            .query_position = 12,
            .pex_leaf_index = 2
        },
        search::seed {
            .sequence = query_span.subspan(18,6),
            .num_errors = 2,
            .query_position = 18,
            .pex_leaf_index = 3
        },
        search::seed {
            .sequence = query_span.subspan(24,6),
            .num_errors = 2,
            .query_position = 24,
            .pex_leaf_index = 4
        }
    };
    EXPECT_EQ(seeds, expected_seeds);
}
