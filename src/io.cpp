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

reference_input read_reference(std::filesystem::path const& reference_sequence_path) {
    auto reference_reader = ivio::fasta::reader{{ .input = reference_sequence_path }};

    std::vector<std::vector<uint8_t>> reference_sequences{};
    std::vector<std::string> reference_tags{};

    for (auto const record_view : reference_reader) {
        std::vector<uint8_t> sequence(record_view.seq.size());
        ivs::convert_char_to_rank<ivs::d_dna4>(record_view.seq, sequence);

        auto const result = ivs::verify_rank(sequence);
        if (result.has_value()) {
            size_t const position = result.value();
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

        reference_sequences.emplace_back(std::move(sequence));
        reference_tags.emplace_back(record_view.id);
    }

    return reference_input { reference_sequences, reference_tags };
}

std::vector<query> read_queries(std::filesystem::path const& queries_path) {
    auto query_reader = ivio::fastq::reader{{ .input = queries_path }};
    std::vector<query> queries{};

    for (auto const record_view : query_reader) {
        std::string const tag(record_view.id);

        std::vector<uint8_t> sequence(record_view.seq.size());
        ivs::convert_char_to_rank<ivs::d_dna4>(record_view.seq, sequence);

        auto const result = ivs::verify_rank(sequence);
        if (result.has_value()) {
            size_t const position = result.value();
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
        }
        
        queries.emplace_back(tag, std::move(sequence));
    }

    return queries;
}

void save_index_and_data(fmindex_with_metadata const& _index, std::filesystem::path const& _index_path) {
    auto ofs     = std::ofstream(_index_path, std::ios::binary);
    auto archive = cereal::BinaryOutputArchive{ofs};
    archive(_index);
}

fmindex_with_metadata load_index_and_data(std::filesystem::path const& _index_path) {
    auto ifs     = std::ifstream(_index_path, std::ios::binary);
    auto archive = cereal::BinaryInputArchive{ifs};
    auto index = fmindex_with_metadata{};
    archive(index);
    return index;
}

} // namespace io
