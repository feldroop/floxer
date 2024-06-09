#include <alignment.hpp>

#include <seqan3/io/sam_file/detail/cigar.hpp>

#include <gtest/gtest.h>

TEST(alignment, small_wrapped_seqan3) {
    using namespace alignment;

    std::vector<uint8_t> reference{ 0, 0, 1, 2, 1, 3, 0, 2, 2, 3, 0, 1 };
    std::vector<uint8_t> query{ 1, 2, 1, 3, 1, 2, 2 };

    alignment_config const config{
        .reference_span_offset = 0,
        .num_allowed_errors = 2,
        .orientation = query_orientation::forward,
        .mode = alignment_mode::verify_and_return_alignment
    };

    auto const result = align(reference, query, config);

    EXPECT_EQ(result.outcome, alignment_outcome::alignment_exists);

    EXPECT_TRUE(result.alignment.has_value());

    EXPECT_EQ(result.alignment.value().num_errors, 1);
    EXPECT_EQ(result.alignment.value().orientation, query_orientation::forward);
    EXPECT_EQ(result.alignment.value().start_in_reference, 2);
    EXPECT_EQ(result.alignment.value().cigar, seqan3::detail::parse_cigar("4=1X2="));
}
