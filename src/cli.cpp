#include <cli.hpp>

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

    options parse_options(int argc, char ** argv) {
        options opt{};
        sharg::parser cli_parser = create_cli_parser(argc, argv, opt);

        try {
            cli_parser.parse();
        }
        catch (sharg::parser_error const & e) {
            std::cerr << "[CLI PARSER ERROR]\n" << e.what() << '\n';
            exit(-1);
        }

        return opt;
    }
}
