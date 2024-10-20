#include <about_floxer.hpp>
#include <floxer_cli.hpp>

#include <cmath>
#include <stdexcept>

#include <spdlog/fmt/fmt.h>
#include <sharg/all.hpp>

namespace cli {

std::filesystem::path const& command_line_input::reference_path() const {
    return reference_path_.value;
}

std::filesystem::path const& command_line_input::queries_path() const {
    return queries_path_.value;
}

std::filesystem::path const& command_line_input::output_path() const {
    return output_path_.value;
}

std::optional<std::filesystem::path> command_line_input::index_path() const {
    if (index_path_.value.empty()) {
        return std::nullopt;
    } else {
        return index_path_.value;
    }
}

std::optional<std::filesystem::path> command_line_input::logfile_path() const {
        if (logfile_path_.value.empty()) {
        return std::nullopt;
    } else {
        return logfile_path_.value;
    }
}

std::optional<size_t> command_line_input::query_num_errors() const {
    if (query_num_errors_.value != std::numeric_limits<size_t>::max()) {
        return query_num_errors_.value;
    } else {
        return std::nullopt;
    }
}

std::optional<double> command_line_input::query_error_probability() const {
    if (std::isnan(query_error_probability_.value)) {
        return std::nullopt;
    } else {
        return query_error_probability_.value;
    }
}

size_t command_line_input::pex_seed_num_errors() const {
    return pex_seed_num_errors_.value;
}

size_t command_line_input::max_num_anchors() const {
    return max_num_anchors_.value;
}

std::string command_line_input::anchor_group_order() const {
    return anchor_group_order_.value;
}

bool command_line_input::bottom_up_pex_tree_building() const {
    return bottom_up_pex_tree_building_.value;
}

bool command_line_input::use_interval_optimization() const {
    return use_interval_optimization_.value;
}

double command_line_input::extra_verification_ratio() const {
    return extra_verification_ratio_.value;
}

double command_line_input::allowed_interval_overlap_ratio() const {
    return allowed_interval_overlap_ratio_.value;
}

bool command_line_input::direct_full_verification() const {
    return direct_full_verification_.value;
}

size_t command_line_input::num_threads() const {
    return num_threads_.value;
}

std::optional<size_t> command_line_input::timeout_seconds() const {
    if (timeout_seconds_.value == 0) {
        return std::nullopt;
    } else {
        return timeout_seconds_.value;
    }
}

std::optional<std::string> command_line_input::stats_target() const {
    return stats_target_.value;
}

std::string command_line_input::command_line_call() const {
    std::vector<std::string> individual_calls{
        "floxer",

        reference_path_.command_line_call(),
        queries_path_.command_line_call(),
        index_path().has_value() ? index_path_.command_line_call() : "",
        output_path_.command_line_call(),
        logfile_path().has_value() ? logfile_path_.command_line_call() : "",

        query_num_errors().has_value() ? query_num_errors_.command_line_call() : "",
        query_error_probability().has_value() ? query_error_probability_.command_line_call() : "",
        pex_seed_num_errors_.command_line_call(),

        max_num_anchors_.command_line_call(),
        anchor_group_order_.command_line_call(),

        bottom_up_pex_tree_building() ? bottom_up_pex_tree_building_.command_line_call() : "",
        use_interval_optimization() ? use_interval_optimization_.command_line_call() : "",
        extra_verification_ratio_.command_line_call(),
        allowed_interval_overlap_ratio_.command_line_call(),
        direct_full_verification() ? direct_full_verification_.command_line_call() : "",

        num_threads_.command_line_call(),
        timeout_seconds().has_value() ? timeout_seconds_.command_line_call() : "",
        stats_target().has_value() ? stats_target_.command_line_call() : ""
    };

    return fmt::format("{}", fmt::join(individual_calls, ""));
}

void command_line_input::validate() const {
    if (!query_num_errors().has_value() && !query_error_probability().has_value()) {
        throw std::runtime_error(
            "Either a fixed number of errors in the query or an error probability must be given."
        );
    }

    if (
        query_num_errors().has_value() &&
        query_num_errors().value() < pex_seed_num_errors()
    ) {
        throw std::runtime_error(
            fmt::format(
                "The number of errors per query ({}) must be greater or equal than the number of errors "
                "in the PEX tree leaves ({}).",
                query_num_errors().value(),
                pex_seed_num_errors()
            )
        );
    }

    if (allowed_interval_overlap_ratio() == 1.0 && !use_interval_optimization()) {
        throw std::runtime_error(
            "You cannot set the allowed interval overlap ratio without activating the interval optimization."
        );
    }
}

void command_line_input::parse_and_validate(int argc, char ** argv) {
    sharg::parser parser{ about_floxer::program_name, argc, argv, sharg::update_notifications::off };

    parser.info.author = about_floxer::author;
    parser.info.description = { about_floxer::long_description };
    parser.info.email = about_floxer::email;
    parser.info.url = about_floxer::url;
    parser.info.short_description = about_floxer::short_description;
    parser.info.synopsis = {
        "./floxer --reference hg38.fasta --query reads.fastq --error-probability 0.07 --output mapped_reads.bam",
    };
    parser.info.version = about_floxer::version;
    parser.info.date = about_floxer::version_date;

    parser.add_option(reference_path_.value, sharg::config{
        .short_id = reference_path_.short_id,
        .long_id = reference_path_.long_id,
        .description = "The reference sequences in which floxer will search the queries, i.e. the haystack."
            "Only valid DNA sequences using [AaCcGgTt] characters are allowed.",
        .required = true,
        .validator = sharg::input_file_validator{
            {
                "fa", "fasta", "fna", "ffn", "fas", "faa", "mpfa", "frn",
                "fa.gz", "fasta.gz", "fna.gz", "ffn.gz", "fas.gz", "faa.gz", "mpfa.gz", "frn.gz"
            }
        }
    });

    parser.add_option(queries_path_.value, sharg::config{
        .short_id = queries_path_.short_id,
        .long_id = queries_path_.long_id,
        .description = "The queries which floxer will search in the reference, i.e. the needles."
            "Queries that contain character other than [AaCcGgTt] are skipped.",
        .required = true,
        .validator = sharg::input_file_validator{{"fq", "fastq", "fq.gz", "fastq.gz"}}
    });

    parser.add_option(index_path_.value, sharg::config{
        .short_id = index_path_.short_id,
        .long_id = index_path_.long_id,
        .description = "The file where the constructed FM-index will be stored for later use. "
            "If the file already exists, the index will be read "
            "from it instead of newly constructed.",
        .default_message = "no index file"
    });

    parser.add_option(output_path_.value, sharg::config{
        .short_id = output_path_.short_id,
        .long_id = output_path_.long_id,
        .description = "The file where the alignment results will be stored.",
        .required = true,
        .validator = sharg::output_file_validator{ sharg::output_file_open_options::open_or_create, {"bam", "sam"}}
    });

    parser.add_option(logfile_path_.value, sharg::config{
        .short_id = logfile_path_.short_id,
        .long_id = logfile_path_.long_id,
        .description = "If a logfile path is given, a rotating logfile will be created "
            "and debug information will be written to it.",
        .default_message = "no logfile"
    });

    parser.add_option(query_num_errors_.value, sharg::config{
        .short_id = query_num_errors_.short_id,
        .long_id = query_num_errors_.long_id,
        .description = "The number of errors allowed in each query. This is only used if no error "
            "probability is given. Either this or an error probability must be given.",
        .default_message = "no default",
        .validator = sharg::arithmetic_range_validator{0, 4096}
    });

    parser.add_option(query_error_probability_.value, sharg::config{
        .short_id = query_error_probability_.short_id,
        .long_id = query_error_probability_.long_id,
        .description = "The error probability in the queries, per base. If this is given, it is used "
            "rather than the fixed number of errors. Either this or a fixed number of errors must be "
            "given.",
        .default_message = "no default",
        .validator = sharg::arithmetic_range_validator{0.00001, 0.99999}
    });

    parser.add_option(pex_seed_num_errors_.value, sharg::config{
        .short_id = pex_seed_num_errors_.short_id,
        .long_id = pex_seed_num_errors_.long_id,
        .description = "The number of errors in the leaves of the PEX tree that are used as seeds. "
            "The sequences will be searched with this parameter using the FM-index.",
        .advanced = true,
        .validator = sharg::arithmetic_range_validator{0, 3}
    });

    parser.add_option(max_num_anchors_.value, sharg::config{
        .short_id = max_num_anchors_.short_id,
        .long_id = max_num_anchors_.long_id,
        .description = "The maximum number of anchors that are located by the FM index. "
            "This should be increased only carfully, as it is detremental for performance and "
            "might not improve accuracy.",
        .advanced = true
    });

    parser.add_option(anchor_group_order_.value, sharg::config{
        .short_id = anchor_group_order_.short_id,
        .long_id = anchor_group_order_.long_id,
        .description = "The way in which anchor groups returned from the FM Index search are ordered. "
            "The first anchor groups in the ordering are more likely to be included for verification.",
        .advanced = true,
        .validator = sharg::value_list_validator{ std::vector{ "errors_first", "count_first", "hybrid" } }
    });

    parser.add_flag(bottom_up_pex_tree_building_.value, sharg::config{
        .short_id = bottom_up_pex_tree_building_.short_id,
        .long_id = bottom_up_pex_tree_building_.long_id,
        .description = "Build PEX trees using a new bottom up strategy.",
        .advanced = true
    });

    parser.add_flag(use_interval_optimization_.value, sharg::config{
        .short_id = use_interval_optimization_.short_id,
        .long_id = use_interval_optimization_.long_id,
        .description = "Keep track of already verified intervals to avoid repeating alignment.",
        .advanced = true
    });

    parser.add_option(extra_verification_ratio_.value, sharg::config{
        .short_id = extra_verification_ratio_.short_id,
        .long_id = extra_verification_ratio_.long_id,
        .description = "How much additional sequence should be verified at the ends of the verification intervals. "
            "This parameter describes the ratio between the original verification intervals "
            "and the additional sequence. Larger values prevent the repeated verification "
            "of mostly overlapping intervals arising from slightly shifted anchors.",
        .advanced = true
    });

    parser.add_option(allowed_interval_overlap_ratio_.value, sharg::config{
        .short_id = allowed_interval_overlap_ratio_.short_id,
        .long_id = allowed_interval_overlap_ratio_.long_id,
        .description = "In the interval optimization, if a previously verified interval overlaps to this relative amount "
            "with a new verification interval, the new interval will not be verified.",
        .advanced = true,
        .validator = sharg::arithmetic_range_validator{0.0000001, 1.0}
    });

    parser.add_flag(direct_full_verification_.value, sharg::config{
        .short_id = direct_full_verification_.short_id,
        .long_id = direct_full_verification_.long_id,
        .description = "Instead of PEX hierarchical verification, directly verify the whole query for every anchor.",
        .advanced = true
    });

    parser.add_option(num_threads_.value, sharg::config{
        .short_id = num_threads_.short_id,
        .long_id = num_threads_.long_id,
        .description = "The number of threads to use in the different steps of the program.",
        .validator = sharg::arithmetic_range_validator{1, 4096}
    });

    parser.add_option(timeout_seconds_.value, sharg::config{
        .short_id = timeout_seconds_.short_id,
        .long_id = timeout_seconds_.long_id,
        .description = "If given, no new alignments will be started after this amount "
            "of seconds and the program will shut down once the already running "
            "alignment jobs have been completed. The index building and input reading does "
            "not count into this.",
        .default_message = "no default",
        .advanced = true
    });

    parser.add_option(stats_target_.value, sharg::config{
        .short_id = stats_target_.short_id,
        .long_id = stats_target_.long_id,
        .description = "Can be the value 'terminal', then number of stats about input, seeding and alignments "
            "will be written to stderr. Otherwise it can be a path to a file and the stats are written to this "
            "location in TOML format.",
        .advanced = true
    });

    parser.parse();

    validate();
}

} // namespace cli
