#include <alignment.hpp>
#include <cli.hpp>
#include <fmindex.hpp>
#include <input.hpp>
#include <output.hpp>
#include <pex.hpp>
#include <search.hpp>

#include <atomic>
#include <chrono>
#include <exception>
#include <filesystem>
#include <locale>
#include <memory>
#include <ranges>
#include <span>
#include <vector>

#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/std.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

#include <ivsigma/ivsigma.h>

void initialize_logger(std::optional<std::filesystem::path> const logfile_path) {
    std::vector<spdlog::sink_ptr> sinks;
    
    auto console_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
    console_sink->set_level(spdlog::level::info);
    std::string const log_pattern = "[thread %t] %+";
    console_sink->set_pattern(log_pattern);

    sinks.push_back(console_sink);

    if (logfile_path.has_value()) {
        auto const max_logfile_size = 1024 * 1024 * 5; // 5 MB
        auto const max_num_logfiles = 3;

        auto file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
            logfile_path.value(),
            max_logfile_size,
            max_num_logfiles
        );
        file_sink->set_level(spdlog::level::trace);
        file_sink->set_pattern(log_pattern);

        sinks.push_back(file_sink);
    }

    auto logger = std::make_shared<spdlog::logger>("floxer", begin(sinks), end(sinks));
    logger->set_level(spdlog::level::trace);
    logger->flush_on(spdlog::level::debug);
    
    spdlog::set_default_logger(logger);
}

std::string format_elapsed_time(spdlog::stopwatch const& stopwatch) {
    auto const elapsed_seconds = stopwatch.elapsed();
    if (elapsed_seconds <= std::chrono::seconds(60)) {
        return fmt::format("{:.7} seconds", elapsed_seconds);
    }

    size_t const all_in_seconds = elapsed_seconds.count();
    size_t const seconds = all_in_seconds % 60;

    size_t const all_in_minutes = all_in_seconds / 60;
    size_t const minutes = all_in_minutes % 60;

    size_t const all_in_hours = all_in_minutes / 60;
    size_t const hours = all_in_hours % 24;

    if (hours > 0) {
        return fmt::format("{}:{:02}:{:02} hours", hours, minutes, seconds);
    } else {
        return fmt::format("{:02}:{:02} minutes", minutes, seconds);
    }
}

int main(int argc, char** argv) {
    // for large number output formatting
    std::locale::global(std::locale("en_US.UTF-8"));

    cli::command_line_input cli_input;
    try {
        cli_input.parse_and_validate(argc, argv);
    } catch (std::exception const & e) {
        fmt::print(stderr, "[CLI PARSER ERROR]\n{}\n", e.what());
        return -1;
    }

    try {
        initialize_logger(cli_input.logfile_path());
    } catch (std::exception const& e) {
        fmt::print(
            "[ERROR] An error occured while trying to set up logging. "
            "Trying to continue without logging output.\n {}\n",
            e.what()
        );
    }

    spdlog::info("successfully parsed CLI input ... starting");

    auto const command_line_call = cli_input.command_line_call();
    spdlog::debug("command line call: {}", command_line_call);

    spdlog::info("reading reference sequences from {}", cli_input.reference_path());

    std::vector<input::reference_record> references;
    try {
        references = input::read_references(cli_input.reference_path());
    } catch (std::exception const& e) {
        spdlog::error(
            "An error occured while trying to read the reference from "
            "the file {}.\n{}",
            cli_input.reference_path(),
            e.what()
        );
        return -1;
    }

    if (references.empty()) {
        spdlog::error(
            "The reference file {} is empty, which is not allowed.\n",
            cli_input.reference_path()
        );
        return -1;
    }

    fmindex index;
    if (cli_input.index_path().has_value() && std::filesystem::exists(cli_input.index_path().value())) {
        auto const index_path = cli_input.index_path().value();

        try {
            spdlog::info("loading index from {}", index_path);
            
            index = input::load_index(index_path);
        } catch (std::exception const& e) {
            spdlog::error(
                "An error occured while trying to load the index from "
                "the file {}.\n{}\n",
                index_path,
                e.what()
            );
            return -1;
        }
    } else {
        spdlog::info(
            "building index with {} thread{}",
            cli_input.num_threads(),
            cli_input.num_threads() == 1 ? "" : "s"
        );

        spdlog::stopwatch build_index_stopwatch;    

        size_t constexpr suffix_array_sampling_rate = 16; 
        index = fmindex(
            references | std::views::transform(&input::reference_record::rank_sequence),
            suffix_array_sampling_rate,
            cli_input.num_threads()
        );

        spdlog::info("building index took {}", format_elapsed_time(build_index_stopwatch));

        if (cli_input.index_path().has_value()) {
            auto const index_path = cli_input.index_path().value();

            try {
                spdlog::info("saving index to {}", index_path);

                output::save_index(index, index_path);
            } catch (std::exception const& e) {
                spdlog::warn(
                    "An error occured while trying to write the index to "
                    "the file {}.\nContinueing without saving the index.\n{}\n",
                    index_path,
                    e.what()
                );
            }
        }
    }

    // not available
    // spdlog::info(
    //     "index memory usage: {:L} bytes",
    //     index.memoryUsage()
    // );

    spdlog::info("reading queries from {}", cli_input.queries_path());

    std::vector<input::query_record> fastq_queries;
    try {
        fastq_queries = input::read_queries(cli_input.queries_path());
    } catch (std::exception const& e) {
        spdlog::error(
            "An error occured while trying to read the queries from "
            "the file {}.\n{}\n",
            cli_input.queries_path(),
            e.what()
        );
        return -1;
    }

    if (fastq_queries.empty()) {
        spdlog::warn(
            "The query file {} is empty.\n", cli_input.queries_path()
        );
    }

    auto sam_output = output::sam_output(
        cli_input.output_path(),
        references,
        command_line_call
    );

    search::search_scheme_cache scheme_cache;
    pex_tree_cache tree_cache;

    spdlog::info(
        "aligning {} queries against {} references with {} thread{} "
        "and writing output file to {}",
        fastq_queries.size(),
        references.size(),
        cli_input.num_threads(),
        cli_input.num_threads() == 1 ? "" : "s",
        cli_input.output_path()
    );

    // workaround for handling errors in threads
    std::atomic_bool encountered_error = false;
    std::vector<std::exception_ptr> exceptions{};

    spdlog::stopwatch aligning_stopwatch;  

    #pragma omp parallel for \
        num_threads(cli_input.num_threads()) \
        default(none) \
        private(tree_cache, scheme_cache) \
        shared(fastq_queries, cli_input, references, index, sam_output, exceptions, encountered_error) \
        schedule(static)
    for (size_t i = 0; i < fastq_queries.size(); ++i) {
        if (encountered_error) {
            continue;
        }

        try {
            auto const& fastq_query = fastq_queries[i];
            size_t const query_num_errors = fastq_query.num_errors_from_user_config(cli_input);

            if (fastq_query.rank_sequence.size() <= query_num_errors) {
                spdlog::warn(
                    "Skipping query {}, because its length of {} is smaller or equal to "
                    "the configured number of errors {}.\n",
                    fastq_query.raw_tag,
                    fastq_query.rank_sequence.size(),
                    query_num_errors
                );

                continue;
            }            
            
            if (query_num_errors < cli_input.pex_seed_num_errors()) {
                spdlog::warn(
                    "Skipping query {}, because using the given error rate {}, it has an allowed "
                    "number of errors of {}, which is smaller than the given number of errors "
                    "in PEX tree leaves of {}.\n",
                    fastq_query.raw_tag,
                    // in this case the error probability must have been given
                    cli_input.query_error_probability().value(),
                    query_num_errors,
                    cli_input.pex_seed_num_errors()
                );

                continue;
            }
            
            spdlog::debug("aligning query: {}", fastq_query.raw_tag);

            auto const tree_config = pex_tree_config {
                .total_query_length = fastq_query.rank_sequence.size(),
                .query_num_errors = query_num_errors,
                .leaf_max_num_errors = cli_input.pex_seed_num_errors()
            };
            auto const& tree = tree_cache.get(tree_config);

            auto alignments = alignment::fastq_query_alignments(references.size());

            bool is_reverse_complement = false;
            tree.search(
                references,
                fastq_query.rank_sequence,
                alignments,
                is_reverse_complement,
                scheme_cache,
                index
            );

            spdlog::trace("after forward: {} alignments", alignments.size());

            auto const reverse_complement_fastq_query_rank_sequence = 
                ivs::reverse_complement_rank<ivs::d_dna4>(fastq_query.rank_sequence);
            is_reverse_complement = true;

            tree.search(
                references,
                reverse_complement_fastq_query_rank_sequence,
                alignments,
                is_reverse_complement,
                scheme_cache,
                index
            );

            spdlog::trace("after reverse complement: {} alignments", alignments.size());

            #pragma omp critical
            sam_output.output_for_query(
                fastq_query,
                references,
                alignments
            );
        } catch (...) {
            #pragma omp critical
            exceptions.emplace_back(std::current_exception());

            encountered_error = true;
        }
    }

    for (auto const& e : exceptions) {
        try {
            std::rethrow_exception(e);
        } catch (std::exception const& e) {
            spdlog::error(
                "An error occured while a thread was aligning reads or writing output to "
                "the file {}.\nThe output file is likely incomplete and invalid.\n{}\n",
                cli_input.output_path(),
                e.what()
            );
        } catch (...) {
            spdlog::error("Unknown error occurred\n");
        }
    } 

    if (!exceptions.empty()) {
        return -1;
    }

    spdlog::info("finished aligning successfully in {}", format_elapsed_time(aligning_stopwatch));

    return 0;
}
