#include <alignment.hpp>
#include <cli.hpp>
#include <fmindex.hpp>
#include <input.hpp>
#include <output.hpp>
#include <pex.hpp>
#include <search.hpp>
#include <statistics.hpp>

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
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

int main(int argc, char** argv) {
    cli::command_line_input cli_input;
    try {
        cli_input.parse_and_validate(argc, argv);
    } catch (std::exception const & e) {
        fmt::print(stderr, "[CLI PARSER ERROR]\n{}\n", e.what());
        return -1;
    }

    try {
        output::initialize_logger(cli_input.logfile_path());
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

    input::references references;
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

    if (references.records.empty()) {
        spdlog::error(
            "The reference file {} is empty, which is not allowed.\n",
            cli_input.reference_path()
        );
        return -1;
    }

    spdlog::info("total reference size: {:L}", references.total_sequence_length);

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
            references.records | std::views::transform(&input::reference_record::rank_sequence),
            suffix_array_sampling_rate,
            cli_input.num_threads()
        );

        spdlog::info("building index took {}", output::format_elapsed_time(build_index_stopwatch));

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

    input::queries queries;
    try {
        queries = input::read_queries(cli_input.queries_path());
    } catch (std::exception const& e) {
        spdlog::error(
            "An error occured while trying to read the queries from "
            "the file {}.\n{}\n",
            cli_input.queries_path(),
            e.what()
        );
        return -1;
    }

    spdlog::info("total query size: {:L}", queries.total_sequence_length);

    auto sam_output = output::sam_output(
        cli_input.output_path(),
        references.records,
        command_line_call
    );

    statistics::search_and_alignment_statistics stats{};
    search::search_scheme_cache scheme_cache;
    pex_tree_cache tree_cache;

    spdlog::info(
        "aligning {} queries against {} references with {} thread{} "
        "and writing output file to {}",
        queries.records.size(),
        references.records.size(),
        cli_input.num_threads(),
        cli_input.num_threads() == 1 ? "" : "s",
        cli_input.output_path()
    );

    // workaround for handling errors in threads
    std::atomic_bool encountered_error = false;
    std::vector<std::exception_ptr> exceptions{};

    spdlog::stopwatch aligning_stopwatch;  
    auto progress_bar = output::progress_bar{ .total_num_events = 8263 };
    progress_bar.start();

    // omp declare reduction(name : type : combining_function) initializer(expression)
    #pragma omp declare reduction( \
            statsReduction : \
            statistics::search_and_alignment_statistics : \
            statistics::combine_stats(omp_out, omp_in) \
        ) \
        initializer (omp_priv=omp_orig)

    #pragma omp parallel for \
        num_threads(cli_input.num_threads()) \
        default(none) \
        private(tree_cache, scheme_cache) \
        shared(queries, cli_input, references, index, \
            sam_output, exceptions, encountered_error, progress_bar) \
        schedule(dynamic) \
        reduction(statsReduction:stats)
    for (size_t query_id = 0; query_id < queries.records.size(); ++query_id) {
        if (encountered_error) {
            continue;
        }

        try {
            auto const& query = queries.records[query_id];
            size_t const query_num_errors = query.num_errors_from_user_config(cli_input);

            if (query.rank_sequence.size() <= query_num_errors) {
                spdlog::warn(
                    "Skipping query {}, because its length of {} is smaller or equal to "
                    "the configured number of errors {}.\n",
                    query.raw_tag,
                    query.rank_sequence.size(),
                    query_num_errors
                );

                continue;
            }            
            
            if (query_num_errors < cli_input.pex_seed_num_errors()) {
                spdlog::warn(
                    "Skipping query {}, because using the given error rate {}, it has an allowed "
                    "number of errors of {}, which is smaller than the given number of errors "
                    "in PEX tree leaves of {}.\n",
                    query.raw_tag,
                    // in this case the error probability must have been given
                    cli_input.query_error_probability().value(),
                    query_num_errors,
                    cli_input.pex_seed_num_errors()
                );

                continue;
            }
            
            spdlog::debug("aligning query: {}", query.raw_tag);

            stats.add_query_length(query.rank_sequence.size());

            auto const tree_config = pex_tree_config {
                .total_query_length = query.rank_sequence.size(),
                .query_num_errors = query_num_errors,
                .leaf_max_num_errors = cli_input.pex_seed_num_errors()
            };
            auto const& tree = tree_cache.get(tree_config);

            auto alignments = tree.align_forward_and_reverse_complement(
                references.records,
                query.rank_sequence,
                index,
                scheme_cache,
                stats
            );

            #pragma omp critical
            sam_output.output_for_query(
                query,
                references.records,
                alignments
            );

            #pragma omp critical
            progress_bar.progress(query_id);
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
                "\nAn error occured while a thread was aligning reads or writing output to "
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

    progress_bar.finish();
    spdlog::info("finished aligning successfully in {}", output::format_elapsed_time(aligning_stopwatch));
    
    if (cli_input.print_stats()) {
        stats.print_all_histograms();
    }

    return 0;
}
