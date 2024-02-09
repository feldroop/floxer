#include <io.hpp>

#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>

#include <fmt/core.h>

namespace io {
    input_data read_inputs(
        std::filesystem::path const& reference_sequence_path,
        std::filesystem::path const& queries_path
    ) {
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

        return input_data { reference_sequences, reference_tags, queries };
    }
};
