#include <gtest/gtest.h>

#include <pex.hpp>

#include <span>

TEST(floxer_test, pex_tree_book_test) {
    pex_tree_config config {
        .total_query_length = 12,
        .query_num_errors = 3,
        .leaf_num_errors = 0
    };

    pex_tree tree(config);

    EXPECT_EQ(tree.leaf_query_length(), 3);

    std::vector<uint8_t> const query{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3};
    
    std::vector<std::vector<uint8_t>> const expected_leafs{
        {0, 0, 0},
        {1, 1, 1},
        {2, 2, 2},
        {3, 3, 3}
    };

    auto const actual_leafs = tree.generate_leaf_queries(query);
    std::vector<std::vector<uint8_t>> actual_leaf_vector{};

    for (auto const & actual_leaf : actual_leafs) {
        actual_leaf_vector.emplace_back(actual_leaf.begin(), actual_leaf.end());
    }

    EXPECT_EQ(actual_leaf_vector, expected_leafs);
}

TEST(floxer_test, pex_tree_leaf_error1_test) {
    pex_tree_config config {
        .total_query_length = 12,
        .query_num_errors = 3,
        .leaf_num_errors = 1
    };

    pex_tree tree(config);

    EXPECT_EQ(tree.leaf_query_length(), 6);

    std::vector<uint8_t> const query{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3};
    
    std::vector<std::vector<uint8_t>> const expected_leafs{
        {0, 0, 0, 1, 1, 1},
        {2, 2, 2, 3, 3, 3}
    };

    auto const actual_leafs = tree.generate_leaf_queries(query);
    std::vector<std::vector<uint8_t>> actual_leaf_vector{};

    for (auto const & actual_leaf : actual_leafs) {
        actual_leaf_vector.emplace_back(actual_leaf.begin(), actual_leaf.end());
    }

    EXPECT_EQ(actual_leaf_vector, expected_leafs);
}
