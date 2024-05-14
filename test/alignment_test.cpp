#include <alignment.hpp>

#include <gtest/gtest.h>

TEST(alignment, basic_wfa2) {
    using namespace alignment;
    set_alignment_backend_global(alignment_backend::wfa2);

    std::vector<uint8_t> reference{ 0, 0, 1, 2, 1, 3, 0, 2, 2, 3, 0, 1 };
    std::vector<uint8_t> query{ 1, 2, 1, 3, 1, 2, 2 };

    alignment_config const config{
        .reference_span_offset = 0,
        .num_allowed_errors = 2,
        .orientation = query_orientation::forward,
        .mode = alignment_mode::verify_and_return_alignment
    };

    aligner wfa2_aligner;
    auto const result = wfa2_aligner.align(reference, query, config);

    EXPECT_EQ(result.outcome, alignment_outcome::alignment_exists);

    EXPECT_TRUE(result.alignment.has_value());

    EXPECT_EQ(result.alignment.value().num_errors, 1);
    EXPECT_EQ(result.alignment.value().orientation, query_orientation::forward);
    EXPECT_EQ(result.alignment.value().start_in_reference, 2);
    // EXPECT_EQ(result.alignment.value().cigar, ...);
}
