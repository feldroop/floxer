#include <cli.hpp>
#include <about_floxer.hpp>

#include <stdexcept>
#include <string>

#include <fmt/core.h>
#include <sharg/all.hpp>

namespace cli {

options parse_and_validate_options(int argc, char ** argv) {
    sharg::parser parser{ "floxer", argc, argv, sharg::update_notifications::off };

    parser.info.author = "Felix Leander Droop";
    parser.info.description = { about_floxer::long_description };
    parser.info.email = "felix.droop@fu-berlin.de";
    parser.info.url = about_floxer::url;
    parser.info.short_description = "FM-index longread PEX-based aligner";
    parser.info.synopsis = {
        "./floxer --reference hg38.fasta --query reads.fastq --index hg38.index "
        "--errors 7 --output mapped_reads.bam",
        "./floxer --index hg38.index --query reads.fastq --errors 7 --output mapped_reads.bam",
    };
    parser.info.version = about_floxer::version;
    parser.info.date = about_floxer::version_date;

    options opt{};

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

    char constexpr query_num_erros_short = 'e';
    std::string const query_num_erros_long = "query-errors";
    parser.add_option(opt.query_num_errors, sharg::config{
        .short_id = query_num_erros_short, 
        .long_id = query_num_erros_long, 
        .description = "The number of errors allowed in each query. This is only used if no error "
            "probability is given. Either this or an error probability must be given.",
            .default_message = "no default",
        .validator = sharg::arithmetic_range_validator{0, 4096}
    });

    char constexpr query_error_probability_short = 'p';
    std::string const query_error_probability_long = "error-probability";
    parser.add_option(opt.query_error_probability, sharg::config{
        .short_id = query_error_probability_short, 
        .long_id = query_error_probability_long,
        .description = "The error probability in the queries, per base. If this is given, it is used "
            "rather than the fixed number of errors. Either this or a fixed number of errors must be "
            "given.",
        .default_message = "no default",
        .validator = sharg::arithmetic_range_validator{0.00001, 0.99999}
    });

    parser.add_option(opt.pex_leaf_num_errors, sharg::config{
        .short_id = 'l', 
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

    parser.parse();

    if (!(
            parser.is_option_set(query_num_erros_short) || 
            parser.is_option_set(query_num_erros_long)
        ) &&
        !(
            parser.is_option_set(query_error_probability_short) || 
            parser.is_option_set(query_error_probability_long)
        )) {
        throw std::runtime_error(
            "Either a fixed number of errors in the query or an error probability must be given."
        );
    }

    return opt;
}

} // namespace cli
