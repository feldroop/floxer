#pragma once

#include <fmindex.hpp>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace input {

struct reference_record {
    std::string const raw_tag;
    std::string const sam_format_sanitized_name;
    std::vector<uint8_t> const sequence;
    size_t const sequence_length;
};

struct query_record {
    std::string const raw_tag;
    std::string const sam_format_sanitized_name;
    std::vector<uint8_t> const sequence;
    std::string const quality;
    size_t const sequence_length;
};

std::vector<reference_record> read_references(std::filesystem::path const& reference_sequence_path);

std::vector<query_record> read_queries(std::filesystem::path const& queries_path);

fmindex load_index(std::filesystem::path const& _index_path);

} // namespace input
