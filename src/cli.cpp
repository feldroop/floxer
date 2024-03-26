#include <cli.hpp>
#include <about_floxer.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

#include <fmt/core.h>
#include <sharg/all.hpp>

namespace cli {

static const char reference_path_option_short_id = 'r';
static const std::string reference_path_option_long_id = "reference";

static const char query_path_option_short_id = 'q';
static const std::string query_path_option_long_id = "queries";

static const char index_path_option_short_id = 'i';
static const std::string index_path_option_long_id = "index";

static const char output_path_option_short_id = 'o';
static const std::string output_path_option_long_id = "output";

static const char query_num_errors_option_short_id = 'e';
static const std::string query_num_errors_option_long_id = "query-errors";

static const char query_error_probability_option_short_id = 'p';
static const std::string query_error_probability_option_long_id = "error-probability";

static const char pex_leaf_num_errors_option_short_id = 'l';
static const std::string pex_leaf_num_errors_option_long_id = "pex-leaf-errors";

static const char num_threads_option_short_id = 't';
static const std::string num_threads_option_long_id = "threads";

bool options::query_num_errors_was_set() const {
    return query_num_errors != std::numeric_limits<size_t>::max();
}

bool options::query_error_probability_was_set() const {
    return !std::isnan(query_error_probability);
}

std::string options::command_line_call() const {
    std::string const index_call = fmt::format(" --{} {}", index_path_option_long_id, index_path.filename().c_str());
    std::string const query_num_errors_call =
        fmt::format(" --{} {}", query_num_errors_option_long_id, query_num_errors);
    std::string const query_error_probability_call =
        fmt::format(" --{} {}", query_error_probability_option_long_id, query_error_probability);
    
    return fmt::format(
        "floxer --{} {} --{} {}{} --{} {}{}{} --{} {} --{} {}",
        reference_path_option_long_id, reference_sequence_path.filename().c_str(),
        query_path_option_long_id, queries_path.filename().c_str(),
        index_path.empty() ? "" : index_call,
        output_path_option_long_id, output_path.filename().c_str(),
        query_num_errors_was_set() ? query_num_errors_call : "",
        query_error_probability_was_set() ? query_error_probability_call : "",
        pex_leaf_num_errors_option_long_id, pex_leaf_num_errors,
        num_threads_option_long_id, num_threads
    );
}

void validate_parsed_options(options const& opt) {
    if (!opt.query_num_errors_was_set() && !opt.query_error_probability_was_set()) {
        throw std::runtime_error(
            "Either a fixed number of errors in the query or an error probability must be given."
        );
    }

    if (!opt.query_error_probability_was_set() && opt.query_num_errors < opt.pex_leaf_num_errors) {
        throw std::runtime_error(
            fmt::format(
                "The number of errors per query ({}) must be greater or equal than the number of errors "
                "in the PEX tree leaves ({}).",
                opt.query_num_errors,
                opt.pex_leaf_num_errors
            )
        );
    }
}

options parse_and_validate_options(int argc, char ** argv) {
    sharg::parser parser{ about_floxer::program_name, argc, argv, sharg::update_notifications::off };

    parser.info.author = about_floxer::author;
    parser.info.description = { about_floxer::long_description };
    parser.info.email = about_floxer::email;
    parser.info.url = about_floxer::url;
    parser.info.short_description = about_floxer::short_description;
    parser.info.synopsis = {
        "./floxer --reference hg38.fasta --query reads.fastq --index hg38.index "
        "--error-probability 0.25 --pex-leaf-errors 2 --output mapped_reads.bam --threads 4",
        "./floxer --reference hg38.fasta --query reads.fastq --query-errors 7 --output mapped_reads.bam",
    };
    parser.info.version = about_floxer::version;
    parser.info.date = about_floxer::version_date;

    options opt{};

    parser.add_option(opt.reference_sequence_path, sharg::config{
        .short_id = reference_path_option_short_id, 
        .long_id = reference_path_option_long_id,
        .description = "The reference sequences in which floxer will search the queries, i.e. the haystack."
            "Only valid DNA sequences using [AaCcGgTt] characters are allowed.",
        .required = true,
        .validator = sharg::input_file_validator{
            {"fa", "fasta", "fna", "ffn", "fas", "faa", "mpfa", "frn"}
        }
    });

    parser.add_option(opt.queries_path, sharg::config{
        .short_id = query_path_option_short_id, 
        .long_id = query_path_option_long_id, 
        .description = "The queries which floxer will search in the reference, i.e. the needles."
            "Queries that contain character other than [AaCcGgTt] are skipped.",
        .required = true,
        .validator = sharg::input_file_validator{{"fq", "fastq"}}
    });

    parser.add_option(opt.index_path, sharg::config{
        .short_id = index_path_option_short_id, 
        .long_id = index_path_option_long_id, 
        .description = "The file where the constructed FM-index will be stored for later use. "
            "If the file already exists, the index will be read "
            "from it instead of newly constructed."
    });

    parser.add_option(opt.output_path, sharg::config{
        .short_id = output_path_option_short_id, 
        .long_id = output_path_option_long_id, 
        .description = "The file where the results will be stored.",
        .required = true,
        .validator = sharg::output_file_validator{{"bam", "sam"}}
    });

    parser.add_option(opt.query_num_errors, sharg::config{
        .short_id = query_num_errors_option_short_id, 
        .long_id = query_num_errors_option_long_id, 
        .description = "The number of errors allowed in each query. This is only used if no error "
            "probability is given. Either this or an error probability must be given.",
            .default_message = "no default",
        .validator = sharg::arithmetic_range_validator{0, 4096}
    });

    parser.add_option(opt.query_error_probability, sharg::config{
        .short_id = query_error_probability_option_short_id, 
        .long_id = query_error_probability_option_long_id,
        .description = "The error probability in the queries, per base. If this is given, it is used "
            "rather than the fixed number of errors. Either this or a fixed number of errors must be "
            "given.",
        .default_message = "no default",
        .validator = sharg::arithmetic_range_validator{0.00001, 0.99999}
    });

    parser.add_option(opt.pex_leaf_num_errors, sharg::config{
        .short_id = pex_leaf_num_errors_option_short_id, 
        .long_id = pex_leaf_num_errors_option_long_id, 
        .description = "The number of errors in the leaves of the PEX tree. "
            "The seed sequences will be searched with this parameter using the FM-index.",
        .validator = sharg::arithmetic_range_validator{0, 3}
    });

    parser.add_option(opt.num_threads, sharg::config{
        .short_id = num_threads_option_short_id, 
        .long_id = num_threads_option_long_id, 
        .description = "The number of threads to use in the different steps of the program.",
        .validator = sharg::arithmetic_range_validator{1, 64}
    });

    parser.parse();

    validate_parsed_options(opt);

    return opt;
}

} // namespace cli
