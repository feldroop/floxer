#pragma once

#include <floxer_cli.hpp>
#include <input.hpp>
#include <intervals.hpp>
#include <mutex_wrapper.hpp>
#include <output.hpp>
#include <pex.hpp>
#include <search.hpp>
#include <statistics.hpp>
#include <atomic>
#include <stdexcept>
#include <memory>
#include <variant>

#include <msd/channel.hpp>
#define BS_THREAD_POOL_ENABLE_PRIORITY
#include <BS_thread_pool.hpp>

namespace parallelization {

struct search_task_result {
    input::query_record query;
    size_t query_index;
    pex::pex_tree pex_tree;

    std::vector<search::anchor_package> anchor_packages;
};

enum class spawning_outcome {
    success, input_error, input_exhausted
};

spawning_outcome spawn_search_task(
    input::queries& queries,
    size_t const query_index,
    cli::command_line_input const& cli_input,
    search::searcher const& searcher,
    mutex_guarded<statistics::search_and_alignment_statistics>& global_stats,
    BS::thread_pool& thread_pool,
    msd::channel<std::optional<search_task_result>>& channel,
    std::atomic_bool& threads_should_stop
);

// this data will be shared between all of the verification tasks using a shared pointer
// some of the member are guarded by mutex/atomic, others are read only
// some of the members are owned by all of the verification tasks, some of them are just references to the main thread data
struct shared_verification_data {
    input::query_record const query;
    size_t const query_index;
    input::references const& references;
    pex::pex_tree const pex_tree;
    pex::pex_verification_config const config;
    intervals::verified_intervals_for_all_references verified_intervals_forward;
    intervals::verified_intervals_for_all_references verified_intervals_reverse_complement;
    mutex_guarded<alignment::query_alignments> alignments;
    mutex_guarded<output::alignment_output>& alignment_output;
    std::atomic_size_t num_verification_tasks_remaining;
    mutex_guarded<statistics::search_and_alignment_statistics>& global_stats;
    std::atomic_bool& threads_should_stop;

    shared_verification_data(
        input::query_record const query_,
        size_t const query_index_,
        input::references const& references_,
        pex::pex_tree const pex_tree_,
        cli::command_line_input const& cli_input,
        mutex_guarded<output::alignment_output>& alignment_output_,
        size_t const num_verification_tasks,
        mutex_guarded<statistics::search_and_alignment_statistics>& global_stats,
        std::atomic_bool& threads_should_stop
    );
};

void spawn_verification_task(
    search::anchor_package package,
    std::shared_ptr<shared_verification_data> data,
    BS::thread_pool& thread_pool
);

} // namespace parallelization
