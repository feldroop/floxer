#include <alignment.hpp>
#include <floxer_cli.hpp>
#include <fmindex.hpp>
#include <input.hpp>
#include <intervals.hpp>
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
#include <mutex>
#include <ranges>
#include <span>
#include <thread>
#include <vector>

#define BS_THREAD_POOL_ENABLE_PRIORITY
#include <BS_thread_pool.hpp>

#include <msd/channel.hpp>

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
        // TODO make this configurable via CLI
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

    input::queries queries(cli_input);

    auto const searcher = search::searcher {
        .index = index,
        .num_reference_sequences = references.records.size(),
        .config = search::search_config{
            .max_num_anchors = cli_input.max_num_anchors(),
            .anchor_group_order = search::anchor_group_order_from_string(cli_input.anchor_group_order())
        }
    };

    pex::pex_alignment_config const pex_alignment_config(searcher, cli_input);

    auto alignment_output = output::create_alignment_output(
        cli_input.output_path(),
        references.records
    );
    std::mutex alignment_output_mutex;

    statistics::search_and_alignment_statistics global_stats{};

    bool all_queries_started = false;
    size_t num_queries_started = 0;
    size_t num_queries_finished = 0;
    std::atomic_bool threads_should_stop = false;

    if (cli_input.timeout_seconds().has_value()) {
        std::thread([&cli_input, &threads_should_stop] {
            std::this_thread::sleep_for(std::chrono::seconds(*cli_input.timeout_seconds()));
            threads_should_stop = true;
        }).detach();
    }

    BS::thread_pool thread_pool(cli_input.num_threads());
    msd::channel<parallelization::align_task_result> channel;

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
    while (num_queries_started < cli_input.num_threads()) {
        auto const spawning_outcome = spawn_alignment_task(
            queries,
            thread_pool,
            threads_should_stop,
            cli_input,
            global_stats,
            num_queries_started,
            references,
            pex_alignment_config,
            alignment_output,
            alignment_output_mutex,
            channel
        );

        if (spawning_outcome == parallelization::spawning_outcome::input_error) {
            thread_pool.purge();
            thread_pool.wait();
            return -1;
        }

        if (spawning_outcome == parallelization::spawning_outcome::input_exhausted) {
            all_queries_started = true;
            break;
        }

        assert(spawning_outcome == parallelization::spawning_outcome::success);
        ++num_queries_started;
    }

    if (num_queries_started == 0) {
        spdlog::warn("empty input file {}", cli_input.queries_path());
        return 0;
    }

    for (parallelization::align_task_result const& res : channel) {
        if (res.is_error()) {
            thread_pool.purge();
            thread_pool.wait();

            spdlog::error(
                "An error occured while a thread was aligning reads or writing output. "
                "The output file is likely incomplete and invalid.\n{}\n",
                res.get_exception().what()
            );

            return -1;
        }

        ++num_queries_finished;
        statistics::combine_stats(global_stats, res.get_success_result());

        if (!all_queries_started) {
            auto const spawning_outcome = spawn_alignment_task(
                queries,
                thread_pool,
                threads_should_stop,
                cli_input,
                global_stats,
                num_queries_started,
                references,
                pex_alignment_config,
                alignment_output,
                alignment_output_mutex,
                channel
            );

            if (spawning_outcome == parallelization::spawning_outcome::input_error) {
                thread_pool.purge();
                thread_pool.wait();
                return -1;
            }

            if (spawning_outcome == parallelization::spawning_outcome::input_exhausted) {
                all_queries_started = true;
            } else {
                assert(spawning_outcome == parallelization::spawning_outcome::success);
                ++num_queries_started;
            }
        }

        if (num_queries_started == num_queries_finished) {
            break;
        }
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
        // TODO remove the stats to terminal printing
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
