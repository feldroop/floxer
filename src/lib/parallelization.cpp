#include <alignment.hpp>
#include <parallelization.hpp>
#include <verification.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/fmt/std.h>

#include <limits>

namespace parallelization {

void spawn_search_task(
    mutex_guarded<input::queries>& queries,
    input::references const& references,
    cli::command_line_input const& cli_input,
    search::searcher const& searcher,
    mutex_guarded<output::alignment_output>& alignment_output,
    mutex_guarded<statistics::search_and_alignment_statistics>& global_stats,
    BS::thread_pool& thread_pool,
    std::atomic_bool& threads_should_stop
) {
    thread_pool.detach_task(
        [
            &queries,
            &references,
            &cli_input,
            &searcher,
            &alignment_output,
            &threads_should_stop,
            &global_stats,
            &thread_pool
        ] {
            if (threads_should_stop) {
                return;
            }

            std::optional<size_t> query_internal_id;

            try {
                std::optional<input::query_record> query_opt;

                {
                    auto && [lock, qs] = queries.lock_unique();
                    query_opt = qs.next();
                }

                if (!query_opt.has_value()) {
                    return;
                }

                auto query = *std::move(query_opt);
                query_internal_id = query.internal_id;

                spdlog::debug("searching query {}: {}", query.internal_id, query.id);

                pex::pex_tree_config const pex_tree_config(query.rank_sequence.size(), cli_input);
                pex::pex_tree const pex_tree(pex_tree_config);

                // maybe parallelize of seeds
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

                spdlog::debug("finished searching query {}: {}", query.internal_id, query.id);

                auto shared_verification_data = std::make_shared<parallelization::shared_verification_data>(
                    std::move(query),
                    references,
                    std::move(pex_tree),
                    cli_input,
                    alignment_output,
                    anchor_packages.size(),
                    global_stats,
                    threads_should_stop
                );

                for (auto& package : anchor_packages) {
                    parallelization::spawn_verification_task(
                        std::move(package),
                        shared_verification_data,
                        thread_pool
                    );
                }

                spawn_search_task(
                    queries,
                    references,
                    cli_input,
                    searcher,
                    alignment_output,
                    global_stats,
                    thread_pool,
                    threads_should_stop
                );
            } catch (std::exception const& e) {
                threads_should_stop = true;
                spdlog::error(
                    "An error occurred while this thread was reading and searching the query no. {}.\n"
                    "Shutting down threads. The output file is likely incomplete. Error message:\n{}",
                    query_internal_id.has_value() ? fmt::format("{}", *query_internal_id) : "<unknown>",
                    e.what()
                );
            }
        },
        BS::pr::low
    );
}

shared_verification_data::shared_verification_data(
    input::query_record const query_,
    input::references const& references_,
    pex::pex_tree const pex_tree_,
    cli::command_line_input const& cli_input,
    mutex_guarded<output::alignment_output>& alignment_output_,
    size_t const num_verification_tasks_,
    mutex_guarded<statistics::search_and_alignment_statistics>& global_stats_,
    std::atomic_bool& threads_should_stop_
) : query{query_},
    references{references_},
    pex_tree{pex_tree_},
    config(cli_input),
    verified_intervals_forward(intervals::create_thread_safe_verified_intervals(
        references.records.size(),
        config.use_interval_optimization,
        config.overlap_rate_that_counts_as_contained
    )),
    verified_intervals_reverse_complement(intervals::create_thread_safe_verified_intervals(
        references.records.size(),
        config.use_interval_optimization,
        config.overlap_rate_that_counts_as_contained
    )),
    alignments(references.records.size()),
    alignment_output{alignment_output_},
    num_verification_tasks_remaining(num_verification_tasks_),
    global_stats{global_stats_},
    threads_should_stop{threads_should_stop_}
{}

void spawn_verification_task(
    search::anchor_package package,
    std::shared_ptr<shared_verification_data> data,
    BS::thread_pool& thread_pool
) {
    thread_pool.detach_task(
        [
            package = std::move(package),
            data
        ] {
            if (data->threads_should_stop) {
                return;
            }

            try {
                spdlog::debug("verifiying package {} of query {}: {}", package.package_id, data->query.internal_id, data->query.id);

                statistics::search_and_alignment_statistics local_stats;

                auto& verified_intervals_for_all_references = package.orientation == alignment::query_orientation::forward ?
                    data->verified_intervals_forward :
                    data->verified_intervals_reverse_complement;

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
                        .already_verified_intervals = verified_intervals_for_all_references.at(anchor.reference_id),
                        .extra_verification_ratio = data->config.extra_verification_ratio,
                        .alignments = data->alignments,
                        .stats = local_stats
                    };

                    verifier.verify(data->config.verification_kind);
                }

                spdlog::debug("finished verifiying package {} of query {}: {}", package.package_id, data->query.internal_id, data->query.id);

                // write to output file if I am the last remaining thread
                if (data->num_verification_tasks_remaining.fetch_sub(1) == 1) {
                    // this locking is only necessary for the mutex wrapper (because this is the last verification task oif this query)
                    auto && [alignments_lock, alignments] = data->alignments.lock_unique();

                    local_stats.add_num_alignments(alignments.size());

                    for (size_t reference_id = 0; reference_id < data->references.records.size(); ++reference_id) {
                        for (auto const& alignment : alignments.to_reference(reference_id)) {
                            local_stats.add_alignment_edit_distance(alignment.num_errors);
                        }
                    }

                    spdlog::debug("(package {}) writing alignments for query {}: {}", package.package_id, data->query.internal_id, data->query.id);

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
                    data->query.internal_id, e.what()
                );
            }
        },
        BS::pr::high
    );
}

} // namespace parallelization