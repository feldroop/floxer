#include <input.hpp>

#include <gtest/gtest.h>

TEST(input, extract_record_id) {
    std::string const id = "kcmieo25789377djs28";
    std::string const tag =  id + " metadata";
    EXPECT_EQ(input::internal::extract_record_id(tag), id);
}

TEST(input, chars_to_rank_sequence_simple) {
    std::string const simple_chars = "ACGTacgt";
    std::vector<uint8_t> const expected_rank_sequence{ 1,2,3,4,1,2,3,4 };
    EXPECT_EQ(input::internal::chars_to_rank_sequence(simple_chars), expected_rank_sequence);
}

TEST(input, chars_to_rank_sequence_sentinel) {
    std::string const chars_with_sentinel = "ACGTacgt$";
    std::vector<uint8_t> const expected_rank_sequence{ 1,2,3,4,1,2,3,4,0 };
    EXPECT_EQ(input::internal::chars_to_rank_sequence(chars_with_sentinel), expected_rank_sequence);
}

TEST(input, chars_to_rank_sequence_invald_chars) {
    std::string const chars_with_invalid = "ACGTacgtW3>"; // 'U' becomes 4, just like 'T'. NOt sure if this is good behavior from ivsigma
    std::vector<uint8_t> const expected_rank_sequence{ 1,2,3,4,1,2,3,4,5,5,5 };
    EXPECT_EQ(input::internal::chars_to_rank_sequence(chars_with_invalid), expected_rank_sequence);
}