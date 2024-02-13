#pragma once

#include <filesystem>

namespace cli {
    struct options {
        std::filesystem::path reference_sequence;
        std::filesystem::path queries;
        std::filesystem::path output_path;
        std::filesystem::path index_path;
        size_t query_num_errors;
        size_t pex_leaf_num_errors = 3;
    };

    options parse_and_validate_options(int argc, char ** argv);
}
