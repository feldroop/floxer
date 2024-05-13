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

TEST(intervals, verified_intervals) {
    using namespace intervals;

    interval_test_cases const cases{};

    // in this this, we'll keep the outcome fixed to test the basic logic
    auto const outcome = alignment::alignment_outcome::alignment_exists;
    auto const found_outcome = std::make_optional(alignment::alignment_outcome::alignment_exists);

    verified_intervals ivls(intervals::use_interval_optimization::on);

    EXPECT_EQ(ivls.size(), 0);

    ivls.insert(cases.ivl1, outcome);
    ivls.insert(cases.ivl2, outcome);

    EXPECT_EQ(ivls.size(), 2);

    EXPECT_EQ(ivls.contains(cases.ivl1), found_outcome);
    EXPECT_EQ(ivls.contains(cases.ivl2), found_outcome);

    // [ivl1), [ivl2)
    EXPECT_EQ(ivls.contains(cases.inside_ivl1), found_outcome);
    EXPECT_EQ(ivls.contains(cases.overlapping_below_ivl1), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.containing_ivl1), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.overlapping_below_ivl2), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.overlapping_above_ivl2), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.between_both), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.overlapping_both), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.containing_both), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.below_both), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.above_both), std::nullopt);

    ivls.insert(cases.ivl3, outcome);
    EXPECT_EQ(ivls.size(), 2);

    // [ivl1+3), [ivl2)
    EXPECT_EQ(ivls.contains(cases.inside_ivl1), found_outcome);
    EXPECT_EQ(ivls.contains(cases.overlapping_below_ivl1), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.containing_ivl1), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.overlapping_below_ivl2), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.overlapping_above_ivl2), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.between_both), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.overlapping_both), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.containing_both), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.below_both), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.above_both), std::nullopt);

    ivls.insert(cases.ivl4, outcome);
    EXPECT_EQ(ivls.size(), 1);

    // [ivl1+2+3+4)
    EXPECT_EQ(ivls.contains(cases.inside_ivl1), found_outcome);
    EXPECT_EQ(ivls.contains(cases.overlapping_below_ivl1), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.containing_ivl1), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.overlapping_below_ivl2), found_outcome);
    EXPECT_EQ(ivls.contains(cases.overlapping_above_ivl2), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.between_both), found_outcome);
    EXPECT_EQ(ivls.contains(cases.overlapping_both), found_outcome);
    EXPECT_EQ(ivls.contains(cases.containing_both), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.below_both), std::nullopt);
    EXPECT_EQ(ivls.contains(cases.above_both), std::nullopt);

    ivls.insert(cases.ivl5, outcome);
    EXPECT_EQ(ivls.size(), 1);

    // [ivl5)
    EXPECT_EQ(ivls.contains(cases.inside_ivl1), found_outcome);
    EXPECT_EQ(ivls.contains(cases.overlapping_below_ivl1), found_outcome);
    EXPECT_EQ(ivls.contains(cases.containing_ivl1), found_outcome);
    EXPECT_EQ(ivls.contains(cases.overlapping_below_ivl2), found_outcome);
    EXPECT_EQ(ivls.contains(cases.overlapping_above_ivl2), found_outcome);
    EXPECT_EQ(ivls.contains(cases.between_both), found_outcome);
    EXPECT_EQ(ivls.contains(cases.overlapping_both), found_outcome);
    EXPECT_EQ(ivls.contains(cases.containing_both), found_outcome);
    EXPECT_EQ(ivls.contains(cases.below_both), found_outcome);
    EXPECT_EQ(ivls.contains(cases.above_both), found_outcome);

    // size should not change when insert repeatedly
    ivls.insert(cases.ivl5, outcome);
    EXPECT_EQ(ivls.size(), 1);
}
