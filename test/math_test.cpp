#include <math.hpp>

#include <gtest/gtest.h>

TEST(math, saturate_value_to_int32_max) {
    size_t const num = 42;
    int32_t const num32 = 42;

    EXPECT_EQ(math::saturate_value_to_int32_max(num), num32);

    size_t const big_num = std::numeric_limits<size_t>::max();
    EXPECT_EQ(math::saturate_value_to_int32_max(big_num), std::numeric_limits<int32_t>::max());

}

TEST(math, ceil_div) {
    EXPECT_EQ(math::ceil_div(100, 8), 13);
    EXPECT_EQ(math::ceil_div(100, 5), 20);
}
