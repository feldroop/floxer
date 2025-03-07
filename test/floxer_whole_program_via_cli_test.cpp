#include <alignment.hpp>

#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <seqan3/io/sam_file/detail/cigar.hpp>
#include <seqan3/io/sam_file/input.hpp>

#include <gtest/gtest.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/std.h>

std::tuple<int, std::string> run_floxer_on_test_data(std::string const& additional_params) {
    std::string output_filename = std::tmpnam(nullptr);
    output_filename += ".sam";
    std::string const command = fmt::format(
        "../bin/floxer "
        "--reference ../test_data/reference.fasta "
        "--queries ../test_data/queries.fastq "
        "--output {} "
        "--interval-optimization "
        "--console-debug-logs "
        "{}",
        output_filename,
        additional_params
    );

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    int const exit_code = std::system(command.c_str());

    return std::make_tuple(exit_code, output_filename);
}

void check_floxer_output_file(std::string const& output_filename) {
    seqan3::sam_file_input floxer_output{output_filename};

    std::unordered_set<std::string> mentioned_query_ids{};

    for (auto const& record : floxer_output) {
        mentioned_query_ids.insert(record.id());

        if (record.id() == "query1" || record.id() == "query6") {
            EXPECT_EQ(record.flag(), seqan3::sam_flag::unmapped);
            continue;
        }

        EXPECT_FALSE(static_cast<bool>(record.flag() & seqan3::sam_flag::unmapped));

        using namespace seqan3::literals;

        // some of the queries have a unique matching position, others have a range of possible positions
        if (record.id() == "query2" && static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand)) {
            EXPECT_EQ(record.reference_position().value(), 48);
            EXPECT_EQ(record.tags().get<"NM"_tag>(), 0);
            EXPECT_EQ(record.cigar_sequence(), seqan3::detail::parse_cigar("12="));
        } else if (record.id() == "query2" && !static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand)) {
            EXPECT_EQ(record.reference_position().value(), 11);
            EXPECT_EQ(record.tags().get<"NM"_tag>(), 0);
            EXPECT_EQ(record.cigar_sequence(), seqan3::detail::parse_cigar("12="));
        } else if (record.id() == "query3" && static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand)) {
            EXPECT_TRUE(record.reference_position().value() >= 17);
            EXPECT_TRUE(record.reference_position().value() <= 26);
            EXPECT_EQ(record.tags().get<"NM"_tag>(), 2);
            EXPECT_EQ(record.cigar_sequence(), seqan3::detail::parse_cigar("6=2I4="));
        } else if (record.id() == "query3" && !static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand)) {
            EXPECT_TRUE(record.reference_position().value() >= 36);
            EXPECT_TRUE(record.reference_position().value() <= 44);
            EXPECT_EQ(record.tags().get<"NM"_tag>(), 2);
            EXPECT_EQ(record.cigar_sequence(), seqan3::detail::parse_cigar("4=2I6="));
        } else if (record.id() == "query4" && static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand)) {
            EXPECT_TRUE(record.reference_position().value() >= 7);
            EXPECT_TRUE(record.reference_position().value() <= 61);
            EXPECT_EQ(record.tags().get<"NM"_tag>(), 2);
            EXPECT_EQ(record.cigar_sequence(), seqan3::detail::parse_cigar("2I10="));
        } else if (record.id() == "query4" && !static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand)) {
            EXPECT_TRUE(record.reference_position().value() >= 54);
            EXPECT_TRUE(record.reference_position().value() <= 61);
            EXPECT_EQ(record.tags().get<"NM"_tag>(), 2);
            EXPECT_EQ(record.cigar_sequence(), seqan3::detail::parse_cigar("10=2I"));
        } else if (record.id() == "query5" && static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand)) {
            EXPECT_EQ(record.reference_position().value(), 53);
            EXPECT_EQ(record.tags().get<"NM"_tag>(), 0);
            EXPECT_EQ(record.cigar_sequence(), seqan3::detail::parse_cigar("12="));
        } else if (record.id() == "query5" && !static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand)) {
            EXPECT_EQ(record.reference_position().value(), 6);
            EXPECT_EQ(record.tags().get<"NM"_tag>(), 0);
            EXPECT_EQ(record.cigar_sequence(), seqan3::detail::parse_cigar("12="));
        }
    }

    std::unordered_set<std::string> expected_query_ids{
        "query1", "query2", "query3", "query4", "query5", "query6"
    };

    EXPECT_EQ(expected_query_ids, mentioned_query_ids);
}

void run_floxer_via_cli_and_check_output(size_t const seed_errors, size_t const num_threads) {
    auto const [exit_code, output_filename] = run_floxer_on_test_data(
        fmt::format(
            "--query-errors 2 "
            "--seed-errors {} "
            "--extra-verification-ratio 2 "
            "--threads {} ",
            seed_errors,
            num_threads
        )
    );

    EXPECT_EQ(exit_code, 0);

    std::string const out = testing::internal::GetCapturedStdout();
    std::string const err = testing::internal::GetCapturedStderr();

    if (exit_code != 0 || !out.empty())  {
        fmt::print("FLOXER STDOUT: {}\n", out);
        fmt::print("FLOXER STDERR: {}\n", err);
    }

    // all of the diagnostic output should be in stderr
    EXPECT_TRUE(out.empty());

    check_floxer_output_file(output_filename);

    std::filesystem::remove(output_filename);
}

TEST(floxer, whole_program_via_cli_old_pex) {
    run_floxer_via_cli_and_check_output(0, 1);
}

TEST(floxer, whole_program_via_cli_adjusted_pex) {
    run_floxer_via_cli_and_check_output(1, 1);
}

TEST(floxer, whole_program_via_cli_multithreaded) {
    run_floxer_via_cli_and_check_output(1, 4);
}
