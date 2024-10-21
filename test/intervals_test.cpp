#include <intervals.hpp>

#include <gtest/gtest.h>

struct interval_test_cases {
    // two disjoint intervals
    intervals::half_open_interval const ivl1{ .start = 5, .end = 11 };
    intervals::half_open_interval const ivl2{ .start = 15, .end = 21 };

    // not quite connects ivl1 and ivl2, merges with ivl1
    intervals::half_open_interval const ivl3{ .start = 11, .end = 14 };

    // connects ivl1+ivl3 and ivl2
    intervals::half_open_interval const ivl4{ .start = 14, .end = 15 };

    // large intervals that covers everything
    intervals::half_open_interval const ivl5{ .start = 0, .end = 100 };

    // some intervals with different relationships to ivl1 and ivl2
    intervals::half_open_interval const inside_ivl1{ .start = 6, .end = 10 };
    intervals::half_open_interval const overlapping_below_ivl1{ .start = 3, .end = 7 };
    intervals::half_open_interval const containing_ivl1{ .start = 3, .end = 14 };

    intervals::half_open_interval const overlapping_below_ivl2{ .start = 13, .end = 18 };
    intervals::half_open_interval const overlapping_above_ivl2{ .start = 17, .end = 23 };

    intervals::half_open_interval const between_both{ .start = 11, .end = 15 };
    intervals::half_open_interval const overlapping_both{ .start = 8, .end = 16 };
    intervals::half_open_interval const containing_both{ .start = 3, .end = 30 };

    intervals::half_open_interval const below_both{ .start = 0, .end = 2 };
    intervals::half_open_interval const above_both{ .start = 22, .end = 24 };
};

TEST(intervals, half_open_interval) {
    using namespace intervals;

    interval_test_cases const cases{};

    EXPECT_EQ(cases.ivl1.relationship_with(cases.inside_ivl1), interval_relationship::contains);
    EXPECT_EQ(cases.ivl1.relationship_with(cases.overlapping_below_ivl1), interval_relationship::overlapping_or_touching_above);
    EXPECT_EQ(cases.ivl1.relationship_with(cases.containing_ivl1), interval_relationship::inside);
    EXPECT_EQ(cases.ivl1.relationship_with(cases.overlapping_below_ivl2), interval_relationship::completely_below);
    EXPECT_EQ(cases.ivl1.relationship_with(cases.overlapping_above_ivl2), interval_relationship::completely_below);
    EXPECT_EQ(cases.ivl1.relationship_with(cases.between_both), interval_relationship::overlapping_or_touching_below);
    EXPECT_EQ(cases.ivl1.relationship_with(cases.overlapping_both), interval_relationship::overlapping_or_touching_below);
    EXPECT_EQ(cases.ivl1.relationship_with(cases.containing_both), interval_relationship::inside);
    EXPECT_EQ(cases.ivl1.relationship_with(cases.below_both), interval_relationship::completely_above);
    EXPECT_EQ(cases.ivl1.relationship_with(cases.above_both), interval_relationship::completely_below);
    EXPECT_EQ(cases.ivl1.relationship_with(cases.ivl1), interval_relationship::equal);

    EXPECT_EQ(cases.ivl2.relationship_with(cases.inside_ivl1), interval_relationship::completely_above);
    EXPECT_EQ(cases.ivl2.relationship_with(cases.overlapping_below_ivl1), interval_relationship::completely_above);
    EXPECT_EQ(cases.ivl2.relationship_with(cases.containing_ivl1), interval_relationship::completely_above);
    EXPECT_EQ(cases.ivl2.relationship_with(cases.overlapping_below_ivl2), interval_relationship::overlapping_or_touching_above);
    EXPECT_EQ(cases.ivl2.relationship_with(cases.overlapping_above_ivl2), interval_relationship::overlapping_or_touching_below);
    EXPECT_EQ(cases.ivl2.relationship_with(cases.between_both), interval_relationship::overlapping_or_touching_above);
    EXPECT_EQ(cases.ivl2.relationship_with(cases.overlapping_both), interval_relationship::overlapping_or_touching_above);
    EXPECT_EQ(cases.ivl2.relationship_with(cases.containing_both), interval_relationship::inside);
    EXPECT_EQ(cases.ivl2.relationship_with(cases.below_both), interval_relationship::completely_above);
    EXPECT_EQ(cases.ivl2.relationship_with(cases.above_both), interval_relationship::completely_below);
    EXPECT_EQ(cases.ivl2.relationship_with(cases.ivl2), interval_relationship::equal);
}

TEST(intervals, half_open_intervals_trim) {
    intervals::half_open_interval const base{ .start = 10, .end = 20 };

    auto const trim0 = base.trim_from_both_sides(0);
    auto const trim1 = base.trim_from_both_sides(1);
    auto const trim5 = base.trim_from_both_sides(5);
    auto const trim10 = base.trim_from_both_sides(10);
    auto const trim25 = base.trim_from_both_sides(25);

    intervals::half_open_interval const expected_trim0{ .start = 10, .end = 20 };
    intervals::half_open_interval const expected_trim1{ .start = 11, .end = 19 };
    intervals::half_open_interval const expected_trim5{ .start = 14, .end = 15 };
    intervals::half_open_interval const expected_trim10{ .start = 10, .end = 11 };
    intervals::half_open_interval const expected_trim25 = expected_trim10;

    EXPECT_EQ(trim0, expected_trim0);
    EXPECT_EQ(trim1, expected_trim1);
    EXPECT_EQ(trim5, expected_trim5);
    EXPECT_EQ(trim10, expected_trim10);
    EXPECT_EQ(trim25, expected_trim25);
}

TEST(intervals, verified_intervals) {
    using namespace intervals;

    interval_test_cases const cases{};

    verified_intervals ivls(use_interval_optimization::on, 1.0);

    ivls.insert(cases.ivl1);
    ivls.insert(cases.ivl2);

    EXPECT_TRUE(ivls.contains(cases.ivl1));
    EXPECT_TRUE(ivls.contains(cases.ivl2));

    // [ivl1), [ivl2)
    EXPECT_TRUE(ivls.contains(cases.inside_ivl1));
    EXPECT_FALSE(ivls.contains(cases.overlapping_below_ivl1));
    EXPECT_FALSE(ivls.contains(cases.containing_ivl1));
    EXPECT_FALSE(ivls.contains(cases.overlapping_below_ivl2));
    EXPECT_FALSE(ivls.contains(cases.overlapping_above_ivl2));
    EXPECT_FALSE(ivls.contains(cases.between_both));
    EXPECT_FALSE(ivls.contains(cases.overlapping_both));
    EXPECT_FALSE(ivls.contains(cases.containing_both));
    EXPECT_FALSE(ivls.contains(cases.below_both));
    EXPECT_FALSE(ivls.contains(cases.above_both));

    ivls.insert(cases.ivl3);

    // [ivl1+3), [ivl2)
    EXPECT_TRUE(ivls.contains(cases.inside_ivl1));
    EXPECT_FALSE(ivls.contains(cases.overlapping_below_ivl1));
    EXPECT_FALSE(ivls.contains(cases.containing_ivl1));
    EXPECT_FALSE(ivls.contains(cases.overlapping_below_ivl2));
    EXPECT_FALSE(ivls.contains(cases.overlapping_above_ivl2));
    EXPECT_FALSE(ivls.contains(cases.between_both));
    EXPECT_FALSE(ivls.contains(cases.overlapping_both));
    EXPECT_FALSE(ivls.contains(cases.containing_both));
    EXPECT_FALSE(ivls.contains(cases.below_both));
    EXPECT_FALSE(ivls.contains(cases.above_both));

    ivls.insert(cases.ivl4);

    // [ivl1+2+3+4)
    EXPECT_TRUE(ivls.contains(cases.inside_ivl1));
    EXPECT_FALSE(ivls.contains(cases.overlapping_below_ivl1));
    EXPECT_FALSE(ivls.contains(cases.containing_ivl1));
    EXPECT_TRUE(ivls.contains(cases.overlapping_below_ivl2)); // !
    EXPECT_FALSE(ivls.contains(cases.overlapping_above_ivl2));
    EXPECT_TRUE(ivls.contains(cases.between_both)); // !
    EXPECT_TRUE(ivls.contains(cases.overlapping_both)); // !
    EXPECT_FALSE(ivls.contains(cases.containing_both));
    EXPECT_FALSE(ivls.contains(cases.below_both));
    EXPECT_FALSE(ivls.contains(cases.above_both));

    ivls.insert(cases.ivl5);

    // [ivl5)
    EXPECT_TRUE(ivls.contains(cases.inside_ivl1));
    EXPECT_TRUE(ivls.contains(cases.overlapping_below_ivl1));
    EXPECT_TRUE(ivls.contains(cases.containing_ivl1));
    EXPECT_TRUE(ivls.contains(cases.overlapping_below_ivl2));
    EXPECT_TRUE(ivls.contains(cases.overlapping_above_ivl2));
    EXPECT_TRUE(ivls.contains(cases.between_both));
    EXPECT_TRUE(ivls.contains(cases.overlapping_both));
    EXPECT_TRUE(ivls.contains(cases.containing_both));
    EXPECT_TRUE(ivls.contains(cases.below_both));
    EXPECT_TRUE(ivls.contains(cases.above_both));

    // size should not change when insert repeatedly
    ivls.insert(cases.ivl5);
}

TEST(intervals, verified_intervals_overlapping) {
    using namespace intervals;

    interval_test_cases const cases{};

    verified_intervals ivls(use_interval_optimization::on, 0.5);

    ivls.insert(cases.ivl1);

    EXPECT_TRUE(ivls.contains(cases.inside_ivl1));
    EXPECT_TRUE(ivls.contains(cases.overlapping_below_ivl1));
    EXPECT_TRUE(ivls.contains(cases.containing_ivl1));
    EXPECT_FALSE(ivls.contains(cases.overlapping_below_ivl2));
    EXPECT_FALSE(ivls.contains(cases.overlapping_above_ivl2));
    EXPECT_FALSE(ivls.contains(cases.between_both));
    EXPECT_FALSE(ivls.contains(cases.overlapping_both));
    EXPECT_FALSE(ivls.contains(cases.containing_both));
    EXPECT_FALSE(ivls.contains(cases.below_both));
    EXPECT_FALSE(ivls.contains(cases.above_both));

    ivls.insert(cases.ivl2);

    EXPECT_TRUE(ivls.contains(cases.inside_ivl1));
    EXPECT_TRUE(ivls.contains(cases.overlapping_below_ivl1));
    EXPECT_TRUE(ivls.contains(cases.containing_ivl1));
    EXPECT_TRUE(ivls.contains(cases.overlapping_below_ivl2));
    EXPECT_TRUE(ivls.contains(cases.overlapping_above_ivl2));
    EXPECT_FALSE(ivls.contains(cases.between_both));
    EXPECT_FALSE(ivls.contains(cases.overlapping_both));
    EXPECT_FALSE(ivls.contains(cases.containing_both));
    EXPECT_FALSE(ivls.contains(cases.below_both));
    EXPECT_FALSE(ivls.contains(cases.above_both));

    ivls.insert(cases.ivl3);
    ivls.insert(cases.ivl4);

    EXPECT_TRUE(ivls.contains(cases.inside_ivl1));
    EXPECT_TRUE(ivls.contains(cases.overlapping_below_ivl1));
    EXPECT_TRUE(ivls.contains(cases.containing_ivl1));
    EXPECT_TRUE(ivls.contains(cases.overlapping_below_ivl2));
    EXPECT_TRUE(ivls.contains(cases.overlapping_above_ivl2));
    EXPECT_TRUE(ivls.contains(cases.between_both));
    EXPECT_TRUE(ivls.contains(cases.overlapping_both));
    EXPECT_TRUE(ivls.contains(cases.containing_both));
    EXPECT_FALSE(ivls.contains(cases.below_both));
    EXPECT_FALSE(ivls.contains(cases.above_both));
}
