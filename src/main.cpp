#include <iostream>
#include <filesystem>

#include <ivio/ivio.h>
#include <sharg/all.hpp>

#include "floxer.cpp"

struct cli_options {
    std::filesystem::path reference_genome;
    std::filesystem::path queries;
    std::filesystem::path output_file;
};

sharg::parser create_cli_parser(int argc, char ** argv, cli_options& opt) {
    sharg::parser parser{"floxer", argc, argv};

    parser.info.author = "Felix Leander Droop";
    parser.info.date = "05.02.2024";
    parser.info.description = {
        "An exact longread aligner using FM-index search with optimal search schemes, "
        "the PEX hierarchical verification scheme and a highly parallel pairwise alignment implementation."
    };
    parser.info.email = "felix.droop@fu-berlin.de";
    parser.info.short_description = "FM-index longread PEX-based aligner";
    parser.info.synopsis = {
        "./floxer --reference hg38.fasta --query reads.fastq --output mapped_reads.bam"
    };
    parser.info.version = "0.0.0";

    parser.add_option(opt.reference_genome, sharg::config{
        .short_id = 'r', 
        .long_id = "reference", 
        .description = "The reference genome in which floxer will search the queries, i.e. the haystack.",
        .required = true,
        .validator = sharg::input_file_validator{{"fa", "fasta", ".fna"}}
    });
    parser.add_option(opt.queries, sharg::config{
        .short_id = 'q', 
        .long_id = "query", 
        .description = "The queries which floxer will search in the reference, i.e. the needles.",
        .required = true,
        .validator = sharg::input_file_validator{{"fq", "fastq"}}
    });
    parser.add_option(opt.output_file, sharg::config{
        .short_id = 'o', 
        .long_id = "output", 
        .description = "The file where the results will be stored.",
        .required = true,
        .validator = sharg::output_file_validator{{"bam", "sam"}}
    });

    return parser;
}

int main(int argc, char** argv) {
    cli_options opt{};
    sharg::parser cli_parser = create_cli_parser(argc, argv, opt);

    try {
        cli_parser.parse();
    }
    catch(sharg::parser_error const & e) {
        std::cerr << "[CLI PARSER ERROR]\n" << e.what() << '\n';
    }
    
    std::cout << "---{ welcome to floxer }---\n\n"
        << "- reference path: " << opt.reference_genome.c_str() << '\n'
        << "- query path: " << opt.queries.c_str() << '\n'
        << "- output path: " << opt.output_file.c_str() << "\n\n";

    pex_tree p(12, 3);
    p.debug_print();

    std::cout << "---{ Reading input files..." << std::endl; 
    
    auto reference_reader = ivio::fasta::reader{{ .input = opt.reference_genome }};

    std::string reference_genome{};
    std::string reference_tags{};

    for (auto record_view : reference_reader) {
        reference_tags += reference_tags.empty() ? "" : ":";
        reference_tags += record_view.id;
        reference_genome += record_view.seq;
    }

    // TODO build FM-index using fmindex-collection
    auto query_reader = ivio::fastq::reader{{ .input = opt.queries }};
    auto queries = std::vector(begin(query_reader), end(query_reader));

    std::cout << "                              ...done }---\n"; 

    return 0;
}
