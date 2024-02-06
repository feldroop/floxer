#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <ivio/ivio.h>

namespace io {
    struct input_data {
        std::string reference_genome;
        std::string reference_tags;
        std::vector<ivio::fastq::record> queries;
    };

    input_data read_inputs(
        std::filesystem::path const& reference_genome_path,
        std::filesystem::path const& queries_path
    );
}
