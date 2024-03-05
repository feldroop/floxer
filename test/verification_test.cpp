#include <gtest/gtest.h>

#include <fmindex.hpp>
#include <pex.hpp>
#include <search.hpp>
#include <verification.hpp>

#include <span>

#include <fmt/core.h>
#include <fmt/ranges.h>
TEST(floxer_test, verification_query_occurs) {
    using namespace verification;

    std::vector<uint8_t> const simple_reference{
        0,0,0,0,0,0,0,0,0,0,
        1,1,1,1,1,1,1,1,1,1,
        2,2,2,2,2,2,2,2,2,2,
        3,3,3,3,3,3,3,3,3,3
    };
    std::vector<uint8_t> const noisy_reference{0,1,3,2,1,0,3,1,0,2,1,0,3,2};

    std::vector<uint8_t> const big_query_perfect_match{0,0,0,0,0,0,1,1,1,1,1,1};
    std::vector<uint8_t> const small_query_perfect_match{1,1,1,1,1};
    std::vector<uint8_t> const tiny_query_perfect_match{0};
    
    std::vector<uint8_t> const big_query_2_errors{0,0,3,0,0,0,1,1,1,2,1,1};
    std::vector<uint8_t> const small_query_1_error{1,0,0,1,1,1};

    std::vector<uint8_t> const noisy_query_2_inserts{1,3,2,1,3,0,3,1,1,0,2};
    std::vector<uint8_t> const noisy_query_2_deletes{3,0,2,0,3,2};
    std::vector<uint8_t> const noisy_query_1_of_each_error_type{1,3,2,0,1,0,1,1,0,1,0,3};
    
    for (size_t num_errors = 0; num_errors < 5; ++ num_errors) {
        auto const result_simple_big_perfect = 
            query_occurs(simple_reference, big_query_perfect_match, num_errors);
        ASSERT_TRUE(result_simple_big_perfect.has_value());
        auto const alignment_simple_big_perfect = result_simple_big_perfect.value();
        ASSERT_EQ(alignment_simple_big_perfect.num_errors, 0);
        ASSERT_EQ(alignment_simple_big_perfect.start_in_reference, 4);
        ASSERT_EQ(alignment_simple_big_perfect.end_in_reference, 16);
        ASSERT_EQ(
            alignment_simple_big_perfect.alignment,
            alignment_from_string("MMMMMMMMMMMM")
        );

        auto const result_simple_small_perfect = 
            query_occurs(simple_reference, small_query_perfect_match, num_errors);
        ASSERT_TRUE(result_simple_small_perfect.has_value());
        auto const alignment_simple_small_perfect = result_simple_small_perfect.value();
        ASSERT_EQ(alignment_simple_small_perfect.num_errors, 0);
        ASSERT_EQ(alignment_simple_small_perfect.start_in_reference, 10);
        ASSERT_EQ(alignment_simple_small_perfect.end_in_reference, 15);
        ASSERT_EQ(
            alignment_simple_small_perfect.alignment,
            alignment_from_string("MMMMM")
        );

        auto const result_simple_tiny_perfect = 
            query_occurs(simple_reference, tiny_query_perfect_match, num_errors);
        ASSERT_TRUE(result_simple_tiny_perfect.has_value());
        auto const alignment_simple_tiny_perfect = result_simple_tiny_perfect.value();
        ASSERT_EQ(alignment_simple_tiny_perfect.num_errors, 0);
        ASSERT_EQ(alignment_simple_tiny_perfect.start_in_reference, 0);
        ASSERT_EQ(alignment_simple_tiny_perfect.end_in_reference, 1);
        ASSERT_EQ(
            alignment_simple_tiny_perfect.alignment,
            alignment_from_string("M")
        );

        auto const result_simple_noisy_nope = 
            query_occurs(simple_reference, noisy_query_2_inserts, num_errors);
        ASSERT_TRUE(!result_simple_noisy_nope.has_value());
    }

    auto const result_simple_big1 = 
        query_occurs(simple_reference, big_query_2_errors, 1);
    ASSERT_TRUE(!result_simple_big1.has_value());
    auto const result_simple_big2 = 
        query_occurs(simple_reference, big_query_2_errors, 2);
    ASSERT_TRUE(result_simple_big2.has_value());
    auto const alignment_simple_big2 = result_simple_big2.value();
    ASSERT_EQ(alignment_simple_big2.num_errors, 2);
    ASSERT_EQ(alignment_simple_big2.start_in_reference, 5);
    ASSERT_EQ(alignment_simple_big2.end_in_reference, 15);
    ASSERT_EQ(
        alignment_simple_big2.alignment,
        alignment_from_string("MMIMMMMMMIMM")
    );

    auto const result_simple_small0 = 
        query_occurs(simple_reference, small_query_1_error, 0);
    ASSERT_TRUE(!result_simple_small0.has_value());
    auto const result_simple_small1 = 
        query_occurs(simple_reference, small_query_1_error, 1);
    ASSERT_TRUE(result_simple_small1.has_value());
    auto const alignment_simple_small1 = result_simple_small1.value();
    ASSERT_EQ(alignment_simple_small1.num_errors, 1);
    ASSERT_EQ(alignment_simple_small1.start_in_reference, 8);
    ASSERT_EQ(alignment_simple_small1.end_in_reference, 13);
    ASSERT_EQ(
        alignment_simple_small1.alignment,
        alignment_from_string("IMMMMM")
    );
    auto const noisy_2_inserts_0_allowed = 
        query_occurs(noisy_reference, noisy_query_2_inserts, 0);
    ASSERT_TRUE(!noisy_2_inserts_0_allowed.has_value());
    auto const noisy_2_inserts_2_allowed = 
        query_occurs(noisy_reference, noisy_query_2_inserts, 2);
    ASSERT_TRUE(noisy_2_inserts_2_allowed.has_value());
    auto const alignment_noisy_2_inserts_2_allowed = noisy_2_inserts_2_allowed.value();
    ASSERT_EQ(alignment_noisy_2_inserts_2_allowed.num_errors, 2);
    ASSERT_EQ(alignment_noisy_2_inserts_2_allowed.start_in_reference, 1);
    ASSERT_EQ(alignment_noisy_2_inserts_2_allowed.end_in_reference, 10);
    ASSERT_EQ(
        alignment_noisy_2_inserts_2_allowed.alignment,
        alignment_from_string("MMMMIMMMIMM")
    );

    auto const noisy_2_deletes_1_allowed = 
        query_occurs(noisy_reference, noisy_query_2_deletes, 1);
    ASSERT_TRUE(!noisy_2_deletes_1_allowed.has_value());
    auto const noisy_2_deletes_2_allowed = 
        query_occurs(noisy_reference, noisy_query_2_deletes, 2);
    ASSERT_TRUE(noisy_2_deletes_2_allowed.has_value());
    auto const alignment_noisy_2_deletes_2_allowed = noisy_2_deletes_2_allowed.value();
    ASSERT_EQ(alignment_noisy_2_deletes_2_allowed.num_errors, 2);
    ASSERT_EQ(alignment_noisy_2_deletes_2_allowed.start_in_reference, 8);
    ASSERT_EQ(alignment_noisy_2_deletes_2_allowed.end_in_reference, 14);
    // one of the deletes can actually be modeled as an insertion
    // this is the currently intended behavior
    ASSERT_EQ(
        alignment_noisy_2_deletes_2_allowed.alignment,
        alignment_from_string("IMMDMMM")
    );

    auto const noisy_all_3 = 
        query_occurs(noisy_reference, noisy_query_1_of_each_error_type, 3);
    ASSERT_TRUE(noisy_all_3.has_value());
    auto const alignment_noisy_all_3 = noisy_all_3.value();
    ASSERT_EQ(alignment_noisy_all_3.num_errors, 3);
    ASSERT_EQ(alignment_noisy_all_3.start_in_reference, 1);
    ASSERT_EQ(alignment_noisy_all_3.end_in_reference, 13);
    ASSERT_EQ(
        alignment_noisy_all_3.alignment,
        alignment_from_string("MMMIMMMMMDMMM")
    );
}
