#include <parallelization.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/fmt/std.h>

namespace parallelization {

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
) {
    input::query_record query;

    try {
        auto query_opt = queries.next();

        if (!query_opt.has_value()) {
            return spawning_outcome::input_exhausted;
        }

        query = *std::move(query_opt);
    } catch (std::exception const& e) {
        spdlog::error(
            "An error occured while trying to read the queries from "
            "the file {}.\n{}\n",
            cli_input.queries_path(),
            e.what()
        );
        threads_should_stop = true;
        return spawning_outcome::input_error;
    }

    stats.add_query_length(query.rank_sequence.size());

    pex::pex_tree_config const pex_tree_config(query.rank_sequence.size(), cli_input);

    thread_pool.detach_task(
        [
            query_index,
            query = std::move(query),
            &references,
            pex_tree_config,
            &pex_alignment_config,
            &alignment_output,
            &alignment_output_mutex,
            &channel,
            &threads_should_stop
        ] {
            if (threads_should_stop) {
                return;
            }

            try {
                statistics::search_and_alignment_statistics local_stats{};

                spdlog::debug("({}) aligning query: {}", query_index, query.id);

                pex::pex_tree const pex_tree(pex_tree_config);

                auto alignments = pex_tree.align_forward_and_reverse_complement(
                    references.records,
                    query.rank_sequence,
                    pex_alignment_config,
                    local_stats
                );

                {
                    const std::lock_guard<std::mutex> lock(alignment_output_mutex);
                    output::output_for_query(
                        alignment_output,
                        query,
                        references.records,
                        std::move(alignments)
                    );
                }

                spdlog::debug("({}) finished aligning query: {}", query_index, query.id);

                channel << parallelization::align_task_result(local_stats);
            } catch (std::exception const& e) {
                threads_should_stop = true;
                channel << parallelization::align_task_result(e);
            }
        }
    );

    return spawning_outcome::success;
}

} // namespace parallelization