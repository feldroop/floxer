#pragma once

#include <cmath>
#include <filesystem>

namespace cli {

struct options {
    std::filesystem::path reference_sequence;
    std::filesystem::path queries;
    std::filesystem::path output_path;
    std::filesystem::path index_path;
    size_t query_num_errors = 0;
    double query_error_probability = NAN;
    size_t pex_leaf_num_errors = 2;
    size_t num_threads = 1;
};

options parse_and_validate_options(int argc, char ** argv);

} // namespace cli
