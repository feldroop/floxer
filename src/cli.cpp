#include <cli.hpp>

#include <fmt/core.h>
#include <sharg/all.hpp>

namespace cli {

sharg::parser create_cli_parser(int argc, char ** argv, options& opt) {
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
        "./floxer --reference hg38.fasta --query reads.fastq --index hg38.index "
        "--errors 7 --output mapped_reads.bam",
        "./floxer --index hg38.index --query reads.fastq --errors 7 --output mapped_reads.bam",
    };
    parser.info.version = "0.0.0";

    parser.add_option(opt.reference_sequence, sharg::config{
        .short_id = 'r', 
        .long_id = "reference",
        .description = "The reference sequences in which floxer will search the queries, i.e. the haystack."
            "Only valid DNA sequences using [AaCcGgTt] characters are allowed.",
        .required = true,
        .validator = sharg::input_file_validator{
            {"fa", "fasta", "fna", "ffn", "fas", "faa", "mpfa", "frn"}
        }
    });
    parser.add_option(opt.queries, sharg::config{
        .short_id = 'q', 
        .long_id = "query", 
        .description = "The queries which floxer will search in the reference, i.e. the needles."
            "Queries that contain character other than [AaCcGgTt] are skipped.",
        .required = true,
        .validator = sharg::input_file_validator{{"fq", "fastq"}}
    });
    parser.add_option(opt.output_path, sharg::config{
        .short_id = 'o', 
        .long_id = "output", 
        .description = "The file where the results will be stored.",
        .required = true,
        .validator = sharg::output_file_validator{{"bam", "sam"}}
    });
    parser.add_option(opt.index_path, sharg::config{
        .short_id = 'i', 
        .long_id = "index", 
        .description = "The file where the constructed FM-index will be stored for later use. "
            "If the file already exists, the index will be read "
            "from it instead of newly constructed."
    });
    parser.add_option(opt.query_num_errors, sharg::config{
        .short_id = 'e', 
        .long_id = "errors", 
        .description = "The number of errors allowed in each query.",
        .required = true,
        .validator = sharg::arithmetic_range_validator{0, 4096}
    });
    parser.add_option(opt.pex_leaf_num_errors, sharg::config{
        .short_id = 'p', 
        .long_id = "pex-leaf-errors", 
        .description = "The number of errors in the leaves of the PEX tree. "
            "The seed sequences will be searched with this parameter using the FM-index.",
        .validator = sharg::arithmetic_range_validator{0, 3}
    });
    parser.add_option(opt.num_threads, sharg::config{
        .short_id = 't', 
        .long_id = "threads", 
        .description = "The number of threads to use in the different steps of the program.",
        .validator = sharg::arithmetic_range_validator{1, 64}
    });

    return parser;
}

options parse_and_validate_options(int argc, char ** argv) {
    options opt{};
    sharg::parser cli_parser = create_cli_parser(argc, argv, opt);

    try {
        cli_parser.parse();
    }
    catch (sharg::parser_error const & e) {
        fmt::print(stderr, "[CLI PARSER ERROR]\n{}\n", e.what());
        exit(-1);
    }

    return opt;
}

} // namespace cli
