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
    size_t query_index;
    input::references const& references;
    pex::pex_tree const pex_tree;
    pex::pex_verification_config const config;
    mutex_guarded<output::alignment_output>& alignment_output;
    mutex_guarded<statistics::search_and_alignment_statistics>& global_stats;
    std::atomic_bool& threads_should_stop;
};

// TODO create constructor for shared_verification_data where alignments and num_remaining are part of it
void spawn_verification_task(
    search::anchor_package package,
    std::shared_ptr<intervals::verified_intervals_for_all_references> verified_intervals_for_all_references,
    std::shared_ptr<shared_verification_data> data,
    std::shared_ptr<mutex_guarded<alignment::query_alignments>> alignments_ptr, // can't be part of shared_verification_data, because it can't be move initialized
    std::shared_ptr<std::atomic_size_t> num_verification_tasks_remaining, // can't be part of shared_verification_data, because it can't be move initialized
    BS::thread_pool& thread_pool
);

} // namespace parallelization
