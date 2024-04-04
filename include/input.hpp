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

std::vector<reference_record> read_references(std::filesystem::path const& reference_sequence_path);

std::vector<query_record> read_queries(std::filesystem::path const& queries_path);

fmindex load_index(std::filesystem::path const& _index_path);

} // namespace input
