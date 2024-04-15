#pragma once

#include <cli.hpp>
#include <fmindex.hpp>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace input {

struct reference_record {
    size_t const id;
    std::string const raw_tag;
    std::string const sam_format_sanitized_name;
    std::vector<uint8_t> const rank_sequence;
};

struct query_record {
    size_t const id;
    std::string const raw_tag;
    std::string const sam_format_sanitized_name;
    std::vector<uint8_t> const rank_sequence;
    std::string const quality;

    size_t num_errors_from_user_config(cli::command_line_input const& cli_input) const;
};

struct references {
    std::vector<reference_record> records;
    size_t total_sequence_length;
};

struct queries {
    std::vector<query_record> records;
    size_t total_sequence_length;
};

references read_references(std::filesystem::path const& reference_sequence_path);

queries read_queries(std::filesystem::path const& queries_path);

fmindex load_index(std::filesystem::path const& _index_path);

namespace internal {

std::string sanitize_reference_name_for_sam(std::string const& reference_name);

std::string sanitize_query_name_for_sam(std::string const& query_name);

constexpr std::array<char, 256> degenerate_to_simple_char_conversion_table();

std::string replace_degenerate_chars(std::string_view const& sequence);

} // namespace internal

} // namespace input
