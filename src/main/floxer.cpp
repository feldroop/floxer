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

    output::initialize_logger(cli_input.logfile_path());

    spdlog::info("successfully parsed CLI input ... starting");

    // TODO instead just summarize config
    auto const command_line_call = cli_input.command_line_call();
    spdlog::debug("command line call: {}", command_line_call);

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

    spdlog::info(
        "total reference size: {}",
        output::format_large_numer(references.total_sequence_length)
    );

    fmindex index;
    if (
        cli_input.index_path().has_value() &&
        std::filesystem::exists(cli_input.index_path().value())
    ) {
        auto const index_path = cli_input.index_path().value();

        try {
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

        spdlog::stopwatch const build_index_stopwatch;    

        size_t constexpr suffix_array_sampling_rate = 16; 
        index = fmindex(
            references.records | std::views::transform(&input::reference_record::rank_sequence),
            suffix_array_sampling_rate,
            cli_input.num_threads()
        );

        spdlog::info(
            "building index took {}",
            output::format_elapsed_time(build_index_stopwatch.elapsed())
        );

        if (cli_input.index_path().has_value()) {
            output::save_index(index, cli_input.index_path().value());
        }
    }

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

    spdlog::info(
        "total query size: {}",
        output::format_large_numer(queries.total_sequence_length)
    );

    auto alignment_output = output::create_alignment_output(
        cli_input.output_path(),
        references.records
    );

    search::search_scheme_cache scheme_cache;
    pex::pex_tree_cache tree_cache;
    
    // setup for workaround for handling errors in threads
    std::atomic_bool threads_should_stop = false;
    std::vector<std::exception_ptr> exceptions{};

    statistics::search_and_alignment_statistics stats{};

    spdlog::stopwatch const aligning_stopwatch;  

    spdlog::info(
        "aligning {} queries against {} references with {} thread{} "
        "and writing output file to {}",
        queries.records.size(),
        references.records.size(),
        cli_input.num_threads(),
        cli_input.num_threads() == 1 ? "" : "s",
        cli_input.output_path()
    );

    // arcane syntax for declaring a custom OpenMP reduction operation for the statistics object
    // omp declare reduction(name : type : combining_function) initializer(expression)
    #pragma omp declare reduction ( \
            statsReduction : \
            statistics::search_and_alignment_statistics : \
            statistics::combine_stats(omp_out, omp_in) \
        ) \
        initializer (omp_priv=omp_orig)

    #pragma omp parallel for \
        num_threads(cli_input.num_threads()) \
        default(none) \
        private(tree_cache, scheme_cache) \
        shared(queries, cli_input, references, index, alignment_output, \
            aligning_stopwatch, exceptions, threads_should_stop) \
        schedule(dynamic) \
        reduction(statsReduction:stats)
    for (size_t query_id = 0; query_id < queries.records.size(); ++query_id) {
        if (threads_should_stop) {
            continue;
        }

        if (cli_input.timeout_seconds().has_value()) {
            auto const timeout = std::chrono::seconds(cli_input.timeout_seconds().value());
            if (aligning_stopwatch.elapsed() >= timeout) {
                threads_should_stop = true;
                continue;
            }
        }

        try {
            auto const& query = queries.records[query_id];
            size_t const query_num_errors = query.num_errors_from_user_config(cli_input);
            stats.add_query_length(query.rank_sequence.size());

            // two cases that likely don't occur in practice where the error are configured in a way such that the 
            // alignment algorithm makes no sense and floxer just flags them as unaligned
            if (
                query.rank_sequence.size() <= query_num_errors ||
                query_num_errors < cli_input.pex_seed_num_errors()
            ) {
                spdlog::debug(
                    "({}/{}) skipping query: {} due to bad num_errors configuration", 
                    query_id, queries.records.size(), query.raw_tag
                );

                auto no_alignments = alignment::query_alignments(references.records.size());

                #pragma omp critical
                output::output_for_query(
                    alignment_output,
                    query,
                    references.records,
                    std::move(no_alignments)
                );

                continue;
            }            
            
            spdlog::debug("({}/{}) aligning query: {}", query_id, queries.records.size(), query.raw_tag);

            auto const searcher = search::searcher{
                .index = index,
                .num_reference_sequences = references.records.size(),
                .scheme_cache = scheme_cache,
                .config = search::search_config{
                    .max_num_raw_anchors = cli_input.max_num_raw_anchors()
                }
            };

            auto const tree_config = pex::pex_tree_config {
                .total_query_length = query.rank_sequence.size(),
                .query_num_errors = query_num_errors,
                .leaf_max_num_errors = cli_input.pex_seed_num_errors()
            };

            auto const& tree = tree_cache.get(tree_config);

            auto alignments = tree.align_forward_and_reverse_complement(
                references.records,
                query.rank_sequence,
                searcher,
                stats
            );

            #pragma omp critical
            output::output_for_query(
                alignment_output,
                query,
                references.records,
                std::move(alignments)
            );

            spdlog::debug("({}/{}) finished aligning query: {}", query_id, queries.records.size(), query.raw_tag);
        } catch (...) {
            #pragma omp critical
            exceptions.emplace_back(std::current_exception());

            threads_should_stop = true;
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

    if (threads_should_stop) {
        spdlog::info(
            "Timed out after {}. Aligned {} queries.",
            output::format_elapsed_time(aligning_stopwatch.elapsed()),
            stats.num_queries()
        );
        return 0;
    }

    spdlog::info(
        "finished aligning successfully in {}",
        output::format_elapsed_time(aligning_stopwatch.elapsed())
    );
    
    if (cli_input.print_stats()) {
        for (auto const& formatted_statistic : stats.format_statistics()) {
            spdlog::info("{}", formatted_statistic);
        }
    }

    return 0;
}
