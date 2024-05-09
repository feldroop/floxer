#pragma once

#include <floxer_cli.hpp>
#include <fmindex.hpp>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace input {

struct reference_record {
    std::string const id;
    std::vector<uint8_t> const rank_sequence;
    size_t const internal_id;
};

struct query_record {
    std::string const id;
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

std::string extract_record_id(std::string_view const& reference_name);

std::vector<uint8_t> chars_to_rank_sequence(std::string_view const chars);

} // namespace internal

} // namespace input
