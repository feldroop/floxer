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

    pex::pex_verification_config const pex_verification_config(cli_input);

    mutex_guarded<output::alignment_output> alignment_output(
        cli_input.output_path(),
        references.records
    );

    mutex_guarded<statistics::search_and_alignment_statistics> global_stats;

    bool all_search_tasks_started = false;
    size_t num_search_tasks_started = 0;
    size_t num_search_tasks_finished = 0;
    std::atomic_bool threads_should_stop = false;

    if (cli_input.timeout_seconds().has_value()) {
        std::thread([&cli_input, &threads_should_stop] {
            std::this_thread::sleep_for(std::chrono::seconds(*cli_input.timeout_seconds()));
            threads_should_stop = true;
            spdlog::warn("Timeout happened. Shutting down threads now. The output file might be incomplete.");
        }).detach();
    }

    BS::thread_pool thread_pool(cli_input.num_threads());

    // std::nullopt sent means an error occurred
    msd::channel<std::optional<parallelization::search_task_result>> search_task_result_channel;

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
    while (num_search_tasks_started < 2 * cli_input.num_threads()) { // TODO experiment with other values than 2
        // TODO multiple queries per search task,
        // do io on thread. how to synchronize in case of exhausted input? another atomic bool needed?
        auto const spawning_outcome = parallelization::spawn_search_task(
            queries,
            cli_input,
            searcher,
            global_stats,
            thread_pool,
            search_task_result_channel,
            threads_should_stop
        );

        if (spawning_outcome == parallelization::spawning_outcome::input_error) {
            thread_pool.purge();
            thread_pool.wait();
            return -1;
        }

        if (spawning_outcome == parallelization::spawning_outcome::input_exhausted) {
            all_search_tasks_started = true;
            break;
        }

        assert(spawning_outcome == parallelization::spawning_outcome::success);
        ++num_search_tasks_started;
    }

    if (num_search_tasks_started == 0) {
        spdlog::warn("empty input file {}", cli_input.queries_path());
        return 0;
    }

    // TODO completely kick out channel and detach veri tasks from inside the search tasks
    while (num_search_tasks_started != num_search_tasks_finished) {
        std::optional<parallelization::search_task_result> res_opt;
        search_task_result_channel >> res_opt;

        // std::nullopt means an error must have happened
        if (!res_opt) {
            thread_pool.purge();
            thread_pool.wait();

            return -1;
        }

        ++num_search_tasks_finished;
        auto res = *std::move(res_opt);

        // TODO construct this already in searcher
        auto shared_verification_data = std::make_shared<parallelization::shared_verification_data>(
            std::move(res.query),
            references,
            std::move(res.pex_tree),
            cli_input,
            alignment_output,
            res.anchor_packages.size(),
            global_stats,
            threads_should_stop
        );

        for (auto& package : res.anchor_packages) {
            parallelization::spawn_verification_task(
                std::move(package),
                shared_verification_data,
                thread_pool
            );
        }
        // TODO decouple query index and num search tasks
        if (!all_search_tasks_started) {
            auto const spawning_outcome = parallelization::spawn_search_task(
                queries,
                cli_input,
                searcher,
                global_stats,
                thread_pool,
                search_task_result_channel,
                threads_should_stop
            );

            if (spawning_outcome == parallelization::spawning_outcome::input_error) {
                thread_pool.purge();
                thread_pool.wait();
                return -1;
            }

            if (spawning_outcome == parallelization::spawning_outcome::input_exhausted) {
                all_search_tasks_started = true;
            } else {
                assert(spawning_outcome == parallelization::spawning_outcome::success);
                ++num_search_tasks_started;
            }
        }
    }

    // wait for all remaining verification tasks to complete
    thread_pool.wait();

    if (!threads_should_stop) {
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
