#pragma once

#include <filesystem>

namespace cli {
    struct options {
        std::filesystem::path reference_genome;
        std::filesystem::path queries;
        std::filesystem::path output_file;
    };

    options parse_options(int argc, char ** argv);
}
