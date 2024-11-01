#include <alignment.hpp>
#include <floxer_cli.hpp>
#include <fmindex.hpp>
#include <input.hpp>
#include <intervals.hpp>
#include <output.hpp>
#include <pex.hpp>
#include <search.hpp>
#include <statistics.hpp>

#include <atomic>
#include <chrono>
#include <exception>
#include <filesystem>
#include <fstream>
#include <locale>
#include <memory>
#include <mutex>
#include <ranges>
#include <span>
#include <vector>

#define BS_THREAD_POOL_ENABLE_PRIORITY
#include <BS_thread_pool.hpp>

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
    spdlog::debug("command line call: {}", cli_input.command_line_call());

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

        // This sampling rate is a trade-off for high speed. It leads to and index size of 11G
        // for the human genome, which should be tolerable in most applications
        size_t constexpr suffix_array_sampling_rate = 4;
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
        queries = input::read_queries(cli_input);
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

    auto const no_alignments = alignment::query_alignments(references.records.size());
    for (auto const& query_with_invalid_config : queries.records_with_invalid_config) {
        output::output_for_query(
            alignment_output,
            query_with_invalid_config,
            references.records,
            no_alignments
        );
    }

    auto const searcher = search::searcher {
        .index = index,
        .num_reference_sequences = references.records.size(),
        .config = search::search_config{
            .max_num_anchors = cli_input.max_num_anchors(),
            .anchor_group_order = search::anchor_group_order_from_string(cli_input.anchor_group_order())
        }
    };

    auto const pex_alignment_config = pex::pex_alignment_config {
        .searcher = searcher,
        .use_interval_optimization = cli_input.use_interval_optimization() ?
            intervals::use_interval_optimization::on :
            intervals::use_interval_optimization::off,
        .verification_kind = cli_input.direct_full_verification() ?
            pex::verification_kind_t::direct_full :
            pex::verification_kind_t::hierarchical,
        .extra_verification_ratio = cli_input.extra_verification_ratio(),
        .overlap_rate_that_counts_as_contained = cli_input.allowed_interval_overlap_ratio()
    };

    statistics::search_and_alignment_statistics global_stats{};

    // setup for workaround for handling errors in threads
    std::atomic_bool threads_should_stop = false;
    std::mutex global_stats_mutex;
    std::mutex alignment_output_mutex;


    spdlog::stopwatch const aligning_stopwatch;

    BS::thread_pool pool(cli_input.num_threads());

    spdlog::info(
        "aligning {} queries against {} references with {} thread{} "
        "and writing output file to {}",
        queries.records.size(),
        references.records.size(),
        cli_input.num_threads(),
        cli_input.num_threads() == 1 ? "" : "s",
        cli_input.output_path()
    );

    BS::multi_future<void> futures = pool.submit_sequence(
        0ul,
        queries.records.size(),
        [
            &queries,
            &cli_input,
            &references,
            &pex_alignment_config,
            &alignment_output,
            &alignment_output_mutex,
            &aligning_stopwatch,
            &threads_should_stop,
            &global_stats,
            &global_stats_mutex
        ] (size_t const query_index) {
            // TODO better stopping mechanism using purge
            if (threads_should_stop) {
                return;
            }

            // TODO think abou tbetter timeout mechanism
            if (cli_input.timeout_seconds().has_value()) {
                auto const timeout = std::chrono::seconds(cli_input.timeout_seconds().value());
                if (aligning_stopwatch.elapsed() >= timeout) {
                    threads_should_stop = true;
                    return;
                }
            }

            auto const& query = queries.records[query_index];
            statistics::search_and_alignment_statistics local_stats{};
            local_stats.add_query_length(query.rank_sequence.size());

            spdlog::debug("({}/{}) aligning query: {}", query_index, queries.records.size(), query.id);

            size_t const query_num_errors = input::num_errors_from_user_config(query.rank_sequence.size(), cli_input);
            auto const pex_tree_config = pex::pex_tree_config {
                .total_query_length = query.rank_sequence.size(),
                .query_num_errors = query_num_errors,
                .leaf_max_num_errors = cli_input.pex_seed_num_errors(),
                .build_strategy = cli_input.bottom_up_pex_tree_building() ?
                    pex::pex_tree_build_strategy::bottom_up :
                    pex::pex_tree_build_strategy::recursive
            };

            pex::pex_tree const pex_tree(pex_tree_config);

            auto alignments = pex_tree.align_forward_and_reverse_complement(
                references.records,
                query.rank_sequence,
                pex_alignment_config,
                local_stats
            );

            if (cli_input.stats_target().has_value()) {
                const std::lock_guard<std::mutex> lock(global_stats_mutex);
                statistics::combine_stats(global_stats, local_stats);
            }

            {
                const std::lock_guard<std::mutex> lock(alignment_output_mutex);
                output::output_for_query(
                    alignment_output,
                    query,
                    references.records,
                    std::move(alignments)
                );
            }

            spdlog::debug("({}/{}) finished aligning query: {}", query_index, queries.records.size(), query.id);
        },
        BS::pr::highest
    );

    try {
        futures.get();
    } catch (std::exception const& e) {
        spdlog::error(
            "An error occured while a thread was aligning reads or writing output to "
            "the file {}.\nThe output file is likely incomplete and invalid.\n{}\n",
            cli_input.output_path(),
            e.what()
        );
        pool.purge();

        return -1;
    } catch (...) {
        spdlog::error("Unknown error occurred in an aligning thread\n");
        pool.purge();

        return -1;
    }

    if (threads_should_stop) {
        spdlog::info(
            "Timed out after {}. Aligned {} queries.",
            output::format_elapsed_time(aligning_stopwatch.elapsed()),
            global_stats.num_queries()
        );
    } else {
        spdlog::info(
            "finished aligning successfully in {}",
            output::format_elapsed_time(aligning_stopwatch.elapsed())
        );
    }

    if (cli_input.stats_target().has_value()) {
        const std::lock_guard<std::mutex> lock(global_stats_mutex); // technically not necessary (all tasks should be done here), but better safe than sorry

        if (cli_input.stats_target().value() == "terminal") {
            for (auto const& formatted_statistic : global_stats.format_statistics_for_stdout()) {
                spdlog::info("{}", formatted_statistic);
            }
        } else {
            std::ofstream out(cli_input.stats_target().value());
            out << global_stats.format_statistics_as_toml();
        }
    }

    return 0;
}
