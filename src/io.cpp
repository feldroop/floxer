#include <io.hpp>

// obsoluete after multi reference sequence refactor
#include <span>

#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>

#include <fmt/core.h>

namespace io {
    input_data read_inputs(
        std::filesystem::path const& reference_sequence_path,
        std::filesystem::path const& queries_path
    ) {
        auto reference_reader = ivio::fasta::reader{{ .input = reference_sequence_path }};

        std::vector<uint8_t> reference_sequence{};
        std::string reference_combined_tags{};

        for (auto const record_view : reference_reader) {
            reference_combined_tags += reference_combined_tags.empty() ? "" : ":";
            reference_combined_tags += record_view.id;

            size_t const begin = reference_sequence.size();
            size_t const length = record_view.seq.size();

            reference_sequence.resize(reference_sequence.size() + length);
            auto const span = std::span(reference_sequence);

            ivs::convert_char_to_rank<ivs::d_dna4>(record_view.seq, span.subspan(begin, length));
        }

        auto const result = ivs::verify_rank(reference_sequence);
        // TODO add the actual invalid character after multi reference sequence refactor
        if (result.has_value()) {
            fmt::println(stderr, "[INPUT ERROR]\nThe reference genome contains invalid characters.");
            exit(-1);
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
                    tag,
                    sequence[position],
                    position
                );
                continue;
            }
            
            queries.emplace_back(tag, std::move(sequence));
        }

        return input_data { reference_sequence, reference_combined_tags, queries };
    }
};
