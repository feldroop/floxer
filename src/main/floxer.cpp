#include <alignment.hpp>
#include <floxer_cli.hpp>
#include <fmindex.hpp>
#include <input.hpp>
#include <intervals.hpp>
#include <mutex_wrapper.hpp>
#include <output.hpp>
#include <parallelization.hpp>
#include <pex.hpp>
#include <search.hpp>
#include <statistics.hpp>

#include <atomic>
#include <cassert>
#include <chrono>
#include <exception>
#include <filesystem>
#include <fstream>
#include <locale>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <thread>
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

    output::initialize_logger(cli_input.logfile_path(), cli_input.console_debug_logs());

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

    fmindex index;
    if (
        cli_input.index_path().has_value() &&
        std::filesystem::exists(cli_input.index_path().value())
    ) {
        auto const index_path = cli_input.index_path().value();
        spdlog::info("loading index from {}", index_path);

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

    mutex_guarded<input::queries> queries(cli_input);

    auto const searcher = search::searcher {
        .index = index,
        .num_reference_sequences = references.records.size(),
        .config = search::search_config{
            .max_num_anchors_hard = cli_input.max_num_anchors_hard(),
            .max_num_anchors_soft = cli_input.max_num_anchors_soft(),
            .anchor_group_order = search::anchor_group_order_from_string(cli_input.anchor_group_order()),
            .anchor_choice_strategy = search::anchor_choice_strategy_from_string(
                cli_input.anchor_choice_strategy()
            ),
            .erase_useless_anchors = !cli_input.dont_erase_useless_anchors()
        }
    };

    mutex_guarded<output::alignment_output> alignment_output(
        cli_input.output_path(),
        references.records
    );

    mutex_guarded<statistics::search_and_alignment_statistics> global_stats(cli_input.stats_input_hint());
    std::atomic_bool threads_should_stop = false;

    if (cli_input.timeout_seconds().has_value()) {
        std::thread([&cli_input, &threads_should_stop] {
            std::this_thread::sleep_for(std::chrono::seconds(*cli_input.timeout_seconds()));
            threads_should_stop = true;
            spdlog::warn("Timeout happened. Shutting down threads now. The output file might be incomplete.");
        }).detach();
    }

    BS::thread_pool thread_pool(cli_input.num_threads());

    auto const query_file_size_bytes = std::filesystem::file_size(cli_input.queries_path());
    spdlog::info(
        "aligning queries from a {} bytes large file against {} references with {} thread{} "
        "and writing output file to {}",
        output::format_large_number(query_file_size_bytes),
        references.records.size(),
        cli_input.num_threads(),
        cli_input.num_threads() == 1 ? "" : "s",
        cli_input.output_path()
    );

    spdlog::stopwatch const aligning_stopwatch;

    // initialize thread pool task queue with a search task for every thread
    for (size_t t = 0; t < cli_input.num_threads(); ++t) {
        parallelization::spawn_search_task(
            queries,
            references,
            cli_input,
            searcher,
            alignment_output,
            global_stats,
            thread_pool,
            threads_should_stop
        );
    }

    // wait for all tasks to complete
    thread_pool.wait();

    if (threads_should_stop) {
        return -1;
    } else{
        spdlog::info(
            "finished aligning successfully in {}",
            output::format_elapsed_time(aligning_stopwatch.elapsed())
        );
    }

    if (cli_input.stats_target().has_value()) {
        auto && [lock, stats] = global_stats.lock_unique();
        if (cli_input.stats_target().value() == "terminal") {
            for (auto const& formatted_statistic : stats.format_statistics_for_stdout()) {
                spdlog::info("{}", formatted_statistic);
            }
        } else {
            std::ofstream out(cli_input.stats_target().value());
            out << stats.format_statistics_as_toml();
        }
    }

    return 0;
}
