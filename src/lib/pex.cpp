#include <math.hpp>
#include <pex.hpp>

#include <cassert>
#include <ranges>

#include <ivsigma/ivsigma.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/spdlog.h>

namespace pex {

size_t pex_tree::node::length_of_query_span() const {
    return query_index_to - query_index_from + 1;
}

bool pex_tree::node::is_root() const {
    return parent_id == null_id;
}

size_t pex_tree::num_leaves() const {
    return leaves.size();
}

pex_tree::node const& pex_tree::root() const {
    auto const& root = inner_nodes.empty() ? leaves.at(0) : inner_nodes.at(0);
    assert(root.is_root());
    return root;
}

pex_tree const& pex_tree_cache::get(pex_tree_config const config) {
    auto [iter, _] = trees.try_emplace(config.total_query_length, config);
    return iter->second;
}

// ------------------------------ PEX tree building + seeding ------------------------------

pex_tree::pex_tree(pex_tree_config const config) 
    : no_error_leaf_query_length{config.total_query_length / (config.query_num_errors + 1)},
    leaf_max_num_errors{config.leaf_max_num_errors} {
    // use 1 based indices until final computation to make sure to match pseudocode
    add_nodes(
        1,
        config.total_query_length,
        config.query_num_errors,
        null_id
    );
}

void pex_tree::add_nodes(
    size_t const query_index_from,
    size_t const query_index_to,
    size_t const num_errors, 
    size_t const parent_id
) {
    // not sure that this name is the correct meaning of this value from the book
    size_t const num_leafs_left = math::ceil_div(num_errors + 1, 2);

    node const curr_node = {
        parent_id,
        query_index_from - 1, // transform to 0-based index
        query_index_to - 1, // transform to 0-based index
        num_errors
    };

    if (num_errors <= leaf_max_num_errors) {
        leaves.emplace_back(curr_node);
    } else {
        size_t const curr_node_id = inner_nodes.size();
        inner_nodes.emplace_back(curr_node);
        size_t const parent_query_span_length = inner_nodes.back().length_of_query_span();

        // the comment below is the implementation according to the book
        // as of now, I don't know why the split is done like that. It leads to
        // an uneven splitting with one super large seed at the end.
        // (maybe there is a bug like a rounding error, 
        // but I did some of the calculations by hand and it seems correct)
        // size_t const query_split_index = query_index_from + num_leafs_left * no_error_leaf_query_length;

        // simply splitting in half is also not a good option, because the errors for both children are not necessarily the same

        size_t const num_errors_for_left_child = (num_leafs_left * num_errors) / (num_errors + 1);
        size_t const num_errors_for_right_child = ((num_errors + 1 - num_leafs_left) * num_errors) / (num_errors + 1);
        
        size_t const num_errors_for_children = num_errors_for_left_child + num_errors_for_right_child;
        
        double fraction_of_leaf_errors_for_left_child;
        switch (num_errors_for_children) {
            case 0:
                fraction_of_leaf_errors_for_left_child = 0.5;
                break;

            case 1:
                fraction_of_leaf_errors_for_left_child = num_errors_for_left_child == 0 ? 1.0 / 3.0 : 2.0 / 3.0;
                break;
            
            default:
                fraction_of_leaf_errors_for_left_child = num_errors_for_left_child / static_cast<double>(num_errors_for_children);
        }

        size_t const query_split_index = query_index_from + static_cast<size_t>(std::round(
            parent_query_span_length * fraction_of_leaf_errors_for_left_child
        ));
        
        add_nodes(
            query_index_from,
            query_split_index - 1,
            num_errors_for_left_child,
            curr_node_id
        );
        add_nodes(
            query_split_index,
            query_index_to,
            num_errors_for_right_child,
            curr_node_id
        );
    }
}

std::vector<search::seed> pex_tree::generate_seeds(
    std::span<const uint8_t> const query
) const {
    std::vector<search::seed> seeds{};
    seeds.reserve(leaves.size());

    for (auto const& leaf : leaves) {
        auto const seed_span = query.subspan(leaf.query_index_from, leaf.length_of_query_span());
        seeds.emplace_back(seed_span, leaf.num_errors);
    }

    return seeds;
}

// ------------------------------ DOT export ------------------------------

std::string pex_tree::node::dot_statement(size_t const id) const {
    return fmt::format(
        "{} [label=\"errors: {}\\nlength: {}\\nrange: [{},{}]\"];\n",
        id,
        num_errors,
        length_of_query_span(),
        query_index_from,
        query_index_to
    );
}

std::string pex_tree::dot_statement() const {
    std::string dot = fmt::format(
        "graph {{\n"
        "label = \"PEX tree for query length {}, {} errors and leaf threshold {} ({} leaves)\";\n"
        "labelloc = \"t\";\n"
        "node [shape=record];\n",
        inner_nodes[0].query_index_to + 1,
        inner_nodes[0].num_errors,
        leaf_max_num_errors,
        num_leaves()
    );
    
    size_t id = 0;
    for (auto const& inner_node : inner_nodes) {
        add_node_to_dot_statement(inner_node, id, dot);
        ++id;
    }
    for (auto const& leaf_node : leaves) {
        add_node_to_dot_statement(leaf_node, id, dot);
        ++id;
    }

    dot += "}\n";

    return dot;
}

void pex_tree::add_node_to_dot_statement(node const& curr_node, size_t const id, std::string& dot) const {
    dot += curr_node.dot_statement(id);
    if (!curr_node.is_root()) {
        dot += fmt::format("{} -- {};\n", id, curr_node.parent_id);
    }
}

// ------------------------------ hierarchical verification + alignment ------------------------------

alignment::query_alignments pex_tree::align_forward_and_reverse_complement(
    std::vector<input::reference_record> const& references,
    std::span<const uint8_t> const query,
    search::searcher const& searcher,
    statistics::search_and_alignment_statistics& stats
) const {
    auto alignments = alignment::query_alignments(references.size());

    align_query_in_given_orientation(
        references,
        query,
        alignments,
        alignment::query_orientation::forward,
        searcher,
        stats
    );

    auto const reverse_complement_query = 
        ivs::reverse_complement_rank<ivs::d_dna4>(query);

    align_query_in_given_orientation(
        references,
        reverse_complement_query,
        alignments,
        alignment::query_orientation::reverse_complement,
        searcher,
        stats
    );

    stats.add_num_alignments(alignments.size());

    for (size_t reference_id = 0; reference_id < references.size(); ++reference_id) {
        for (auto const& alignment : alignments.to_reference(reference_id)) {
            stats.add_alignment_edit_distance(alignment.num_errors);
        }
    }

    return alignments;
}

void pex_tree::align_query_in_given_orientation(
    std::vector<input::reference_record> const& references,
    std::span<const uint8_t> const query,
    alignment::query_alignments& alignments,
    alignment::query_orientation const orientation,
    search::searcher const& searcher,
    statistics::search_and_alignment_statistics& stats
) const {
    auto const seeds = generate_seeds(query);
    stats.add_statistics_for_seeds(seeds);

    auto const search_result = searcher.search_seeds(seeds);
    stats.add_statistics_for_search_result(search_result);

    for (size_t seed_id = 0; seed_id < seeds.size(); ++seed_id) {
        auto const& anchors_of_seed = search_result.anchors_by_seed[seed_id];
        
        if (anchors_of_seed.excluded) {
            continue;
        }
        
        for (size_t reference_id = 0; reference_id < references.size(); ++reference_id) {
            intervals::interval_set already_verified_intervals;

            for (auto const& anchor : anchors_of_seed.anchors_by_reference[reference_id]) {
                hierarchical_verification(
                    anchor,
                    seed_id,
                    query,
                    orientation,
                    references[reference_id],
                    already_verified_intervals,
                    alignments,
                    stats
                );
            }
        }
    }
}

void pex_tree::hierarchical_verification(
    search::anchor const& anchor,
    size_t const seed_id,
    std::span<const uint8_t> const query,
    alignment::query_orientation const orientation,
    input::reference_record const& reference,
    intervals::interval_set& already_verified_intervals,
    alignment::query_alignments& alignments,
    statistics::search_and_alignment_statistics& stats
) const {
    // this depends on the implementation of generate_leave_queries returning the
    // leaf queries in the same order as the leaves (which it should always do!)
    auto pex_node = leaves.at(seed_id);
    size_t const seed_query_index_from = pex_node.query_index_from;

    size_t const no_extra_wiggle_room = 0;
    auto const root_reference_span_config = internal::compute_reference_span_start_and_length(
        anchor,
        root(),
        seed_query_index_from,
        reference.rank_sequence.size(),
        no_extra_wiggle_room
    );

    if (already_verified_intervals.contains(root_reference_span_config.as_half_open_interval())) {
        // we have already verified the interval where the whole query could be found, according to this anchor
        stats.add_reference_span_size_avoided_root(root_reference_span_config.length);
        return;
    }

    // case for when the whole PEX tree is just a single root
    if (pex_node.is_root()) {
        [[maybe_unused]] auto const outcome = internal::try_to_align_pex_node_query_with_reference_span(
            pex_node,
            reference,
            root_reference_span_config,
            query,
            orientation,
            alignments,
            stats
        );
        assert(outcome == alignment::alignment_outcome::alignment_exists);

        already_verified_intervals.insert(root_reference_span_config.as_half_open_interval());

        return;
    }

    pex_node = inner_nodes.at(pex_node.parent_id);

    while (true) {
        size_t const extra_wiggle_room = 5; // TODO: test different values here
        auto const reference_span_config = internal::compute_reference_span_start_and_length(
            anchor,
            pex_node,
            seed_query_index_from,
            reference.rank_sequence.size(),
            extra_wiggle_room
        );

        auto const outcome = internal::try_to_align_pex_node_query_with_reference_span(
            pex_node,
            reference,
            reference_span_config,
            query,
            orientation,
            alignments,
            stats
        );

        already_verified_intervals.insert(reference_span_config.as_half_open_interval());

        if (outcome == alignment::alignment_outcome::no_adequate_alignment_exists || pex_node.is_root()) {
            break;
        }

        pex_node = inner_nodes.at(pex_node.parent_id);
    }
}

namespace internal {

intervals::half_open_interval span_config::as_half_open_interval() const {
    return intervals::half_open_interval{
        .start = offset,
        .end = offset + length
    };
}

span_config compute_reference_span_start_and_length(
    search::anchor const& anchor,
    pex_tree::node const& pex_node,
    size_t const leaf_query_index_from,
    size_t const full_reference_length,
    size_t const extra_wiggle_room
) {
    int64_t const start_signed = static_cast<int64_t>(anchor.position) - 
        (leaf_query_index_from - pex_node.query_index_from)
        - pex_node.num_errors - extra_wiggle_room;
    size_t const reference_span_start = start_signed >= 0 ? start_signed : 0;
    size_t const reference_span_length = std::min(
        pex_node.length_of_query_span() + 2 * pex_node.num_errors + 1 + 2 * extra_wiggle_room,
        full_reference_length - reference_span_start
    );

    return span_config{
        .offset = reference_span_start,
        .length = reference_span_length
    };
}

alignment::alignment_outcome try_to_align_pex_node_query_with_reference_span(
    pex_tree::node const& pex_node,
    input::reference_record const& reference,
    span_config const reference_span_config,
    std::span<const uint8_t> const query,
    alignment::query_orientation const orientation,
    alignment::query_alignments& alignments,
    statistics::search_and_alignment_statistics& stats
) {
    auto const this_node_query_span = query.subspan(
        pex_node.query_index_from,
        pex_node.length_of_query_span()
    );

    auto const reference_subspan = std::span<const uint8_t>(reference.rank_sequence).subspan(
        reference_span_config.offset,
        reference_span_config.length
    );

    auto const config = alignment::alignment_config {
        .reference_span_offset = reference_span_config.offset,
        .num_allowed_errors = pex_node.num_errors,
        .orientation = orientation,
        .mode = pex_node.is_root() ?
            alignment::alignment_mode::verify_and_return_alignment :
            alignment::alignment_mode::only_verify_existance
    };

    auto const alignment_result = alignment::align(
        reference_subspan,
        this_node_query_span,
        config
    );
    
    if (alignment_result.alignment.has_value()) {
        assert(pex_node.is_root());
        assert(alignment_result.outcome == alignment::alignment_outcome::alignment_exists);

        alignments.insert(
            std::move(alignment_result.alignment.value()),
            reference.id
        );
    }

    if (pex_node.is_root()) {
        stats.add_reference_span_size_aligned_root(reference_span_config.length);
    } else {
        stats.add_reference_span_size_aligned_inner_node(reference_span_config.length);
    }

    return alignment_result.outcome;
}

} // namespace internal

} // namespace pex
