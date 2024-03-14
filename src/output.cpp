#include <output.hpp>

#include <fstream>

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <fmt/core.h>

namespace output {

void save_index(fmindex const& _index, std::filesystem::path const& _index_path) {
    auto ofs     = std::ofstream(_index_path, std::ios::binary);
    auto archive = cereal::BinaryOutputArchive{ofs};
    archive(_index);
}

sam_header::sam_header(std::vector<input::reference_record> const& reference_records) {
    for (auto const& record : reference_records) {
        uint32_t reference_length;
        if (record.sequence_length > std::numeric_limits<uint32_t>::max()) {
            fmt::println(
                "[OUTPUT WARNING]: the sequence {} is too long for the SAM file format (length {})\n"
                "Values in the output file will be set to UINT32_MAX.",
                record.raw_tag,
                record.sequence_length
            );

            reference_length = std::numeric_limits<uint32_t>::max();
        } else {
            reference_length = static_cast<uint32_t>(record.sequence_length);
        }
        
        reference_sequences.push_back(reference_sequence_metadata {
            .name = record.sam_format_sanitized_name,
            .length = reference_length
        });
    }
}

std::string sam_header::to_string() const {
    std::stringstream stream{};

    file_level_metadata.append_to_stream(stream);
    for (auto const& reference_sequence : reference_sequences) {
        reference_sequence.append_to_stream(stream);
    }
    program.append_to_stream(stream);

    return stream.str();
}

} // output
