#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace io {
    struct query {
        std::string const tag;
        std::vector<uint8_t> const sequence;
    };

    struct input_data {
        std::vector<uint8_t> const reference_sequence;
        std::string const reference_combined_tags;
        std::vector<query> const queries;
    };

    input_data read_inputs(
        std::filesystem::path const& reference_sequence_path,
        std::filesystem::path const& queries_path
    );
}
