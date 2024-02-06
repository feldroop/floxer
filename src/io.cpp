#include <filesystem>
#include <string>
#include <vector>

#include <ivio/ivio.h>

class io {
public:
    struct input_data {
        std::string reference_genome;
        std::string reference_tags;
        std::vector<ivio::fastq::record> queries;
    };

    static input_data read_inputs(
        std::filesystem::path const& reference_genome_path,
        std::filesystem::path const& queries_path
    ) {
        auto reference_reader = ivio::fasta::reader{{ .input = reference_genome_path }};

        std::string reference_genome{};
        std::string reference_tags{};

        for (auto record_view : reference_reader) {
            reference_tags += reference_tags.empty() ? "" : ":";
            reference_tags += record_view.id;
            reference_genome += record_view.seq;
        }

        auto query_reader = ivio::fastq::reader{{ .input = queries_path }};
        auto queries = std::vector(begin(query_reader), end(query_reader));

        return input_data { reference_genome, reference_tags, queries };
    }
};
