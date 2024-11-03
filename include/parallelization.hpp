#pragma once

#include <floxer_cli.hpp>
#include <input.hpp>
#include <output.hpp>
#include <pex.hpp>
#include <statistics.hpp>

#include <atomic>
#include <stdexcept>
#include <variant>

#include <msd/channel.hpp>
#define BS_THREAD_POOL_ENABLE_PRIORITY
#include <BS_thread_pool.hpp>

namespace parallelization {

template<typename T>
class task_result {
public:
    task_result() : result(std::runtime_error("empty task result")) {}
    task_result(T success_result) : result(std::move(success_result)) {}
    task_result(std::exception exc) : result(std::move(exc)) {}

    bool is_error() const {
        return std::holds_alternative<std::exception>(result);
    }

    bool is_success() const {
        return !is_error();
    }

    T const& get_success_result() const {
        return std::get<0>(result);
    }
    std::exception const& get_exception() const {
        return std::get<1>(result);
    }

private:
    std::variant<T, std::exception> result;
};

using align_task_result = task_result<statistics::search_and_alignment_statistics>;

enum class spawning_outcome {
    success, input_error, input_exhausted
};

spawning_outcome spawn_alignment_task(
    input::queries& queries,
    BS::thread_pool& thread_pool,
    std::atomic_bool& threads_should_stop,
    cli::command_line_input const& cli_input,
    statistics::search_and_alignment_statistics& stats,
    size_t const query_index,
    input::references const& references,
    pex::pex_alignment_config const& pex_alignment_config,
    output::alignment_output_t& alignment_output,
    std::mutex& alignment_output_mutex,
    msd::channel<parallelization::align_task_result>& channel
);

} // namespace parallelization
