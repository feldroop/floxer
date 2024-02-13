#pragma once

#include <floxer_fmindex.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace io {
    struct query {
        std::string const tag;
        std::vector<uint8_t> const sequence;
    };

    struct reference_input {
        std::vector<std::vector<uint8_t>> const sequences;
        std::vector<std::string> const tags;
    };

    std::vector<query> read_queries(std::filesystem::path const& queries_path);

    reference_input read_reference(std::filesystem::path const& reference_sequence_path);

    void save_index_and_data(FloxerFMIndexWithMetaData const& _index, std::filesystem::path const& _index_path);

    FloxerFMIndexWithMetaData load_index_and_data(std::filesystem::path const& _index_path);
}
