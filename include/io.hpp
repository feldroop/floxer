#pragma once

#include <fmindex.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace io {

struct record {
    std::string const tag;
    std::vector<uint8_t> const sequence;
};

std::vector<record> read_queries(std::filesystem::path const& queries_path);

std::vector<record> read_references(std::filesystem::path const& reference_sequence_path);

void save_index(fmindex const& _index, std::filesystem::path const& _index_path);

fmindex load_index(std::filesystem::path const& _index_path);

} // namespace io
