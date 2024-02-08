#include <io.hpp>

#include <iostream>
#include <span>

#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>

namespace io {
    input_data read_inputs(
        std::filesystem::path const& reference_genome_path,
        std::filesystem::path const& queries_path
    ) {
        auto reference_reader = ivio::fasta::reader{{ .input = reference_genome_path }};

        std::vector<uint8_t> reference_genome{};
        std::string reference_combined_tags{};

        for (auto const record_view : reference_reader) {
            reference_combined_tags += reference_combined_tags.empty() ? "" : ":";
            reference_combined_tags += record_view.id;

            size_t const begin = reference_genome.size();
            size_t const length = record_view.seq.size();

            reference_genome.resize(reference_genome.size() + length);
            auto const span = std::span(reference_genome);

            ivs::convert_char_to_rank<ivs::d_dna4>(record_view.seq, span.subspan(begin, length));
        }

        auto const result = ivs::verify_rank(reference_genome);
        if (result.has_value()) {
            std::cerr << "[INPUT ERROR]\n The reference genome contains invalid characters." << std::endl;
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
                std::cerr << "[INPUT WARNING]\nSkipped the query " << tag
                    << " due to the invalid character " << sequence[position] 
                    << " at position " << position << "." << std::endl;
                continue;
            }
            
            queries.emplace_back(tag, std::move(sequence));
        }

        return input_data { reference_genome, reference_combined_tags, queries };
    }
};
