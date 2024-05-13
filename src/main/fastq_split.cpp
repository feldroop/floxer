#include <about_floxer.hpp>
#include <input.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_set>

#include <seqan3/alphabet/quality/phred94.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <sharg/all.hpp>

std::unordered_set<std::string> read_split_id_file(std::filesystem::path const& split_ids_path) {
    std::ifstream in(split_ids_path);

    std::string line;
    std::unordered_set<std::string> chosen_ids;

    while(std::getline(in, line)) {
        chosen_ids.insert(line);
    }

    return chosen_ids;
}

struct longread_input_traits : seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna15;
    using quality_alphabet = seqan3::phred94;
};

int main(int argc, char** argv) {
    sharg::parser parser{ "fastq_split", argc, argv, sharg::update_notifications::off };

    parser.info.author = about_floxer::author;
    parser.info.description = {
        "Given a list of read IDs and a fastq file, create two new fastq files. "
        "One file containing all the reads with the given IDs from the file and one with the rest."
    };
    parser.info.email = about_floxer::email;
    parser.info.url = about_floxer::url;
    parser.info.short_description = "Split fastq file by given list of IDs";
    parser.info.synopsis = {
        "./fastq_split --split-ids ids.txt --input reads.fastq --chosen-ids-output ids.fastq --rest-output rest.fastq",
    };
    parser.info.version = "1.0.0";
    parser.info.date = "13.05.2024";

    std::filesystem::path split_ids_path{};
    std::filesystem::path input_path{};
    std::filesystem::path chosen_ids_output_path{};
    std::filesystem::path rest_output_path{};

    parser.add_option(split_ids_path, sharg::config{
        .short_id = 's',
        .long_id = "split-ids",
        .description = "The read IDs that should be in one of the two output files.",
        .required = true,
        .validator = sharg::input_file_validator{}
    });

    parser.add_option(input_path, sharg::config{
        .short_id = 'i',
        .long_id = "input",
        .description = "File containing the input reads.",
        .required = true,
        .validator = sharg::input_file_validator{{"fq", "fastq", "fq.gz", "fastq.gz"}}
    });

    parser.add_option(chosen_ids_output_path, sharg::config{
        .short_id = 'c',
        .long_id = "chosen-ids-output",
        .description = "The path where a fastq file with the read IDs from the ID file will be created.",
        .required = true,
        .validator = sharg::output_file_validator{{"fq", "fastq"}}
    });

    parser.add_option(rest_output_path, sharg::config{
        .short_id = 'r',
        .long_id = "rest-output",
        .description = "The path where a fastq file with all of the read IDs "
            "except the ones from the ID file will be created.",
        .required = true,
        .validator = sharg::output_file_validator{{"fq", "fastq"}}
    });

    parser.parse();

    auto const chosen_ids = read_split_id_file(split_ids_path);

    //seqan3::sequen

    seqan3::sequence_file_input<longread_input_traits> input{input_path};
    seqan3::sequence_file_output chosen_ids_output{chosen_ids_output_path};
    seqan3::sequence_file_output rest_output{rest_output_path};

    for (auto const& record : input) {
        auto const id = input::internal::extract_record_id(record.id());
        if (chosen_ids.contains(id)) {
            chosen_ids_output.push_back(record);
        } else {
            rest_output.push_back(record);
        }
    }

    return 0;
}
