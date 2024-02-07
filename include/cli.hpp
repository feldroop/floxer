#pragma once

#include <filesystem>

namespace cli {
    struct options {
        std::filesystem::path reference_genome;
        std::filesystem::path queries;
        std::filesystem::path output_file;
        size_t query_num_errors;
        size_t pex_leaf_num_errors = 3;
    };

    options parse_options(int argc, char ** argv);
}
