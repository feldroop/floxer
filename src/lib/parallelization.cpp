#include <alignment.hpp>
#include <parallelization.hpp>
#include <verification.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/fmt/std.h>

namespace parallelization {

spawning_outcome spawn_search_task(
    input::queries& queries,
    size_t const query_index,
    cli::command_line_input const& cli_input,
    search::searcher const& searcher,
    mutex_guarded<statistics::search_and_alignment_statistics>& global_stats,
    BS::thread_pool& thread_pool,
    msd::channel<std::optional<search_task_result>>& channel,
    std::atomic_bool& threads_should_stop
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

    thread_pool.detach_task(
        [
            &cli_input,
            query = std::move(query),
            query_index,
            &searcher,
            &channel,
            &threads_should_stop,
            &global_stats
        ] {
            if (threads_should_stop) {
                return;
            }

            try {
                spdlog::debug("searching query {}: {}", query_index, query.id);

                pex::pex_tree_config const pex_tree_config(query.rank_sequence.size(), cli_input);
                pex::pex_tree const pex_tree(pex_tree_config);

                auto forward_seeds = pex_tree.generate_seeds(query.rank_sequence);
                auto reverse_complement_seeds = pex_tree.generate_seeds(query.reverse_complement_rank_sequence);

                auto forward_search_result = searcher.search_seeds(forward_seeds);
                auto reverse_complement_search_result = searcher.search_seeds(reverse_complement_seeds);

                std::vector<search::anchor_package> anchor_packages;
                forward_search_result.append_anchor_packages(
                    anchor_packages,
                    cli_input.num_anchors_per_verification_task(),
                    alignment::query_orientation::forward
                );
                reverse_complement_search_result.append_anchor_packages(
                    anchor_packages,
                    cli_input.num_anchors_per_verification_task(),
                    alignment::query_orientation::reverse_complement
                );

                // even if no anchors are found, one empty anchor package is created
                // such that one verification task writes the query as unaligned
                if (anchor_packages.empty()) {
                    anchor_packages.emplace_back(search::anchor_package {
                        .package_id = 0,
                        .anchors{},
                        .orientation = alignment::query_orientation::forward
                    });
                }

                // this is confusing, because the stats are written once for forward and once for reverse complement
                // however the alignment stats are written for everything at once (TODO find a good way to fix this)
                statistics::search_and_alignment_statistics local_stats;
                local_stats.add_query_length(query.rank_sequence.size());
                local_stats.add_statistics_for_seeds(forward_seeds);
                local_stats.add_statistics_for_seeds(reverse_complement_seeds);
                local_stats.add_statistics_for_search_result(forward_search_result);
                local_stats.add_statistics_for_search_result(reverse_complement_search_result);
                {
                    auto && [lock, ref] = global_stats.lock_unique();
                    ref.merge_other_into_this(local_stats);
                }

                spdlog::debug("finished searching query {}: {}", query_index, query.id);

                channel << std::make_optional(search_task_result {
                    .query = std::move(query),
                    .query_index = query_index,
                    .pex_tree = std::move(pex_tree),
                    .anchor_packages = std::move(anchor_packages)
                });
            } catch (std::exception const& e) {
                threads_should_stop = true;
                spdlog::error(
                    "An error occurred while this thread was searching the query no. {}.\n"
                    "Shutting down threads. The output file is likely incomplete. Error message:\n{}",
                    query_index, e.what()
                );
                channel << std::optional<search_task_result>(); // std::nullopt
            }
        },
        BS::pr::low
    );

    return spawning_outcome::success;
}

void spawn_verification_task(
    search::anchor_package package,
    std::shared_ptr<intervals::verified_intervals_for_all_references> verified_intervals_for_all_references,
    std::shared_ptr<shared_verification_data> data,
    std::shared_ptr<mutex_guarded<alignment::query_alignments>> alignments_ptr,
    std::shared_ptr<std::atomic_size_t> num_verification_tasks_remaining,
    BS::thread_pool& thread_pool
) {
    thread_pool.detach_task(
        [
            package = std::move(package),
            verified_intervals_for_all_references,
            alignments_ptr,
            num_verification_tasks_remaining,
            data
        ] {
            if (data->threads_should_stop) {
                return;
            }

            try {
                spdlog::debug("verifiying package {} of query {}: {}", package.package_id, data->query_index, data->query.id);

                statistics::search_and_alignment_statistics local_stats;

                for (auto const anchor : package.anchors) {
                    auto const& pex_leaf_node = data->pex_tree.get_leaves().at(anchor.pex_leaf_index);

                    verification::query_verifier verifier {
                        .pex_tree = data->pex_tree,
                        .anchor = anchor,
                        .pex_leaf_node = pex_leaf_node,
                        .query = package.orientation == alignment::query_orientation::forward ?
                            data->query.rank_sequence : data->query.reverse_complement_rank_sequence,
                        .orientation = package.orientation,
                        .reference = data->references.records[anchor.reference_id],
                        .already_verified_intervals = verified_intervals_for_all_references->at(anchor.reference_id),
                        .extra_verification_ratio = data->config.extra_verification_ratio,
                        .alignments = *alignments_ptr,
                        .stats = local_stats
                    };

                    verifier.verify(data->config.verification_kind);
                }

                spdlog::debug("finished verifiying package {} of query {}: {}", package.package_id, data->query_index, data->query.id);

                // write to output file if I am the last remaining thread
                if (num_verification_tasks_remaining->fetch_sub(1) == 1) {
                    // this locking is only necessary for the mutex wrapper (because this is the last verification task oif this query)
                    auto && [alignments_lock, alignments] = alignments_ptr->lock_unique();

                    local_stats.add_num_alignments(alignments.size());

                    for (size_t reference_id = 0; reference_id < data->references.records.size(); ++reference_id) {
                        for (auto const& alignment : alignments.to_reference(reference_id)) {
                            local_stats.add_alignment_edit_distance(alignment.num_errors);
                        }
                    }

                    spdlog::debug("(package {}) writing alignments for query {}: {}", package.package_id, data->query_index, data->query.id);

                    auto && [output_lock, output] = data->alignment_output.lock_unique();
                    output.write_alignments_for_query(data->query, alignments);
                }

                {
                    auto && [lock, global] = data->global_stats.lock_unique();
                    global.merge_other_into_this(local_stats);
                }
            } catch (std::exception const& e) {
                data->threads_should_stop = true;
                spdlog::error(
                    "An error occurred while this thread was verifying (aligning) the query no. {}.\n"
                    "Shutting down threads. The output file is likely incomplete. Error message:\n{}",
                    data->query_index, e.what()
                );
            }
        },
        BS::pr::high
    );
}

} // namespace parallelization