#include <io.hpp>

#include <fstream>

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>

#include <fmt/core.h>

namespace io {

template<typename Reader>
std::vector<record> read_input_records(Reader && reader, bool const is_queries = false) {
    std::vector<record> records{};

    for (auto const record_view : reader) {
        std::string const tag(record_view.id);
        std::vector<uint8_t> const sequence = ivs::convert_char_to_rank<ivs::d_dna4>(record_view.seq);

        auto const result = ivs::verify_rank(sequence);
        if (result.has_value()) {
            size_t const position = result.value();

            if (is_queries) {
                fmt::print(
                    stderr,
                    "[INPUT WARNING]\nSkipped the query {} "
                    "due to the invalid character {} "
                    "at position {}.\n",
                    record_view.id,
                    sequence[position],
                    position
                );

                continue;
            } else {
                fmt::print(
                    stderr, 
                    "[INPUT ERROR]\nThe reference sequence {} "
                    "contians the invalid character {} "
                    "at position {}.\n",
                    record_view.id,
                    sequence[position],
                    position
                );

                exit(-1);
            }
        }

        records.emplace_back(tag, std::move(sequence));
    }

    return records;
}

std::vector<record> read_references(std::filesystem::path const& reference_sequence_path) {
    return read_input_records(ivio::fasta::reader{{ .input = reference_sequence_path }});
}

std::vector<record> read_queries(std::filesystem::path const& queries_path) {
    return read_input_records(ivio::fastq::reader{{ .input = queries_path }}, true);
}

void save_index(fmindex const& _index, std::filesystem::path const& _index_path) {
    auto ofs     = std::ofstream(_index_path, std::ios::binary);
    auto archive = cereal::BinaryOutputArchive{ofs};
    archive(_index);
}

fmindex load_index(std::filesystem::path const& _index_path) {
    auto ifs     = std::ifstream(_index_path, std::ios::binary);
    auto archive = cereal::BinaryInputArchive{ifs};
    auto index = fmindex{};
    archive(index);
    return index;
}

} // namespace io
