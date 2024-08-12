#include <math.hpp>
#include <pex.hpp>
#include <verification.hpp>

#include <cassert>
#include <ranges>
#include <stdexcept>

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

pex_tree::node const& pex_tree::get_parent_of_child(pex_tree::node const& child) const {
    if (child.is_root()) {
        throw std::runtime_error("tried to get parent of PEX tree root");
    }

    return inner_nodes.at(child.parent_id);
}

pex_tree const& pex_tree_cache::get(pex_tree_config const config) {
    auto [iter, _] = trees.try_emplace(config.total_query_length, config);
    return iter->second;
}

// ------------------------------ PEX tree building + seeding ------------------------------

pex_tree::pex_tree(pex_tree_config const config)
    : no_error_seed_length{config.total_query_length / (config.query_num_errors + 1)},
    leaf_max_num_errors{config.leaf_max_num_errors} {
    switch (config.build_strategy) {
        case pex_tree_build_strategy::recursive:
            // use 1 based indices until final computation to make sure to match pseudocode
            add_nodes_recursive(
                1,
                config.total_query_length,
                config.query_num_errors,
                null_id
            );
            break;
        case pex_tree_build_strategy::bottom_up:
            add_nodes_bottom_up(config);
            break;
        default:
            throw std::runtime_error("(should be unreachable) internal bug in the pex tree construction - build strategy");
    }

    assert(root().num_errors == config.query_num_errors);
    assert(root().query_index_from == 0);
    assert(root().query_index_to == config.total_query_length - 1);
}

void pex_tree::add_nodes_recursive(
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

        // this way of splitting from the book leads to a large seed at the rightmost leaf of the tree.
        // the reason is that the total_query_length is not divisible by (query_num_errors + 1) and the floor
        // of the quotient is chosen for no_error_seed_length. Hence, the whole remainder is covered by
        // the very righmost leaf.
        size_t const query_split_index = query_index_from + num_leafs_left * no_error_seed_length;

        // simply splitting in half is also not a good option, because the errors for both children are not necessarily the same

        size_t const num_errors_for_left_child = (num_leafs_left * num_errors) / (num_errors + 1);
        size_t const num_errors_for_right_child = ((num_errors + 1 - num_leafs_left) * num_errors) / (num_errors + 1);

        add_nodes_recursive(
            query_index_from,
            query_split_index - 1,
            num_errors_for_left_child,
            curr_node_id
        );
        add_nodes_recursive(
            query_split_index,
            query_index_to,
            num_errors_for_right_child,
            curr_node_id
        );
    }
}

void pex_tree::add_nodes_bottom_up(pex_tree_config const& config) {
    // this formula follows from Lemma 1 of the book chapter
    size_t const base_leaf_weight = config.leaf_max_num_errors + 1; // a_i from book
    size_t const num_desired_leaves = math::ceil_div(config.query_num_errors + 1, base_leaf_weight);

    // edge case where tree is only a root
    if (num_desired_leaves == 1) {
        inner_nodes.emplace_back(node {
            .parent_id = null_id,
            .query_index_from = 0,
            .query_index_to = config.total_query_length - 1,
            .num_errors = config.query_num_errors
        });

        return;
    }

    // If we round up in the divison for num_desired_leaves, we actually allow too many errors
    // when setting the number of allowed errors of all leaves to config.leaf_max_num_errors.
    // To avoid this inefficiency, we reduce the number of allowed error in some of the leaves.
    size_t const num_leaves_with_less_errors = ((config.query_num_errors + 1) % base_leaf_weight) > 0 ?
        base_leaf_weight - ((config.query_num_errors + 1) % base_leaf_weight) :
        0;

    create_leaves(config, num_desired_leaves, num_leaves_with_less_errors);

    // this reserve is NECESSARY to avoid the vector reallocating and moving element when growing,
    // because we keep a reference to it (current_level_nodes) while calling emplace_back
    inner_nodes.reserve(num_desired_leaves);
    // secure a position for the root, because it MUST be at index 0
    inner_nodes.emplace_back();

    // merge leaves until the tree is complete
    auto current_level_nodes = std::span(leaves);

    while (current_level_nodes.size() > 3) {
        for (size_t i = 0; i < current_level_nodes.size(); i += 2) {
            size_t const num_remaining_nodes = current_level_nodes.size() - i;

            if (num_remaining_nodes == 1) {
                break;
            }

            // if there is an odd number of nodes on this level, merge the last 3 together
            size_t const num_children = (num_remaining_nodes == 3) ? 3 : 2;

            auto child_nodes = current_level_nodes.subspan(i, num_children);
            size_t const new_parent_id = inner_nodes.size();

            inner_nodes.emplace_back(
                create_parent_node(child_nodes, new_parent_id)
            );
        }

        current_level_nodes = std::span(inner_nodes).last(current_level_nodes.size() / 2);
    }

    inner_nodes.front() = create_parent_node(current_level_nodes, 0);
    inner_nodes.front().parent_id = null_id;

    // free reserved space that is not needed
    inner_nodes.shrink_to_fit();
}

void pex_tree::create_leaves(
    pex_tree_config const& config,
    size_t const num_desired_leaves,
    size_t const num_leaves_with_less_errors
) {
    size_t const base_seed_length = config.total_query_length / num_desired_leaves;
    size_t const seed_length_remainder = config.total_query_length % num_desired_leaves;

    leaves.reserve(num_desired_leaves);

    size_t current_start_index = 0;
    for (size_t i = 0; i < num_desired_leaves; ++i) {
        size_t const curr_leaf_length = i < seed_length_remainder ?
            base_seed_length + 1 :
            base_seed_length;

        size_t const num_errors = i < num_leaves_with_less_errors ?
            config.leaf_max_num_errors - 1 :
            config.leaf_max_num_errors;

        leaves.emplace_back(node{
            .parent_id = 0, // will be set later
            .query_index_from = current_start_index,
            .query_index_to = current_start_index + curr_leaf_length - 1,
            .num_errors = num_errors
        });

        current_start_index += curr_leaf_length;
    }
}

pex_tree::node pex_tree::create_parent_node(std::span<pex_tree::node> const child_nodes, size_t const parent_id) {
    assert(!child_nodes.empty());

    size_t children_errors = 0;
    for (auto& child : child_nodes) {
        child.parent_id = parent_id;
        children_errors += child.num_errors;
    }

    return node {
        .parent_id = 0, // will be set later
        .query_index_from = child_nodes.front().query_index_from,
        .query_index_to = child_nodes.back().query_index_to,
        .num_errors = children_errors + child_nodes.size() - 1
    };
}

std::vector<search::seed> pex_tree::generate_seeds(
    std::span<const uint8_t> const query
) const {
    std::vector<search::seed> seeds{};
    seeds.reserve(leaves.size());

    for (auto const& leaf : leaves) {
        auto const seed_span = query.subspan(leaf.query_index_from, leaf.length_of_query_span());
        seeds.emplace_back(search::seed {
            .sequence = seed_span,
            .num_errors = leaf.num_errors,
            .query_position = leaf.query_index_from
        });
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
        root().query_index_to + 1,
        root().num_errors,
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
    pex_alignment_config const config,
    statistics::search_and_alignment_statistics& stats
) const {
    auto alignments = alignment::query_alignments(references.size());

    align_query_in_given_orientation(
        references,
        query,
        alignments,
        alignment::query_orientation::forward,
        config,
        stats
    );

    auto const reverse_complement_query =
        ivs::reverse_complement_rank<ivs::d_dna4>(query);

    align_query_in_given_orientation(
        references,
        reverse_complement_query,
        alignments,
        alignment::query_orientation::reverse_complement,
        config,
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
    pex_alignment_config const config,
    statistics::search_and_alignment_statistics& stats
) const {
    auto const seeds = generate_seeds(query);
    stats.add_statistics_for_seeds(seeds);

    auto const search_result = config.searcher.search_seeds(seeds);
    stats.add_statistics_for_search_result(search_result);

    std::vector<intervals::verified_intervals> already_verified_intervals_per_reference(
        references.size(), intervals::verified_intervals(config.use_interval_optimization)
    );

    for (size_t seed_id = 0; seed_id < seeds.size(); ++seed_id) {
        auto const& anchors_of_seed = search_result.anchors_by_seed[seed_id];

        if (anchors_of_seed.status == search::seed_status::fully_excluded) {
            continue;
        }

        for (size_t reference_id = 0; reference_id < references.size(); ++reference_id) {
            for (auto const& anchor : anchors_of_seed.anchors_by_reference[reference_id]) {
                // this depends on the implementation of generate_leave_queries returning the
                // leaf queries in the same order as the leaves (which it should always do!)
                auto pex_node = leaves.at(seed_id);

                verification::query_verifier verifier {
                    .pex_tree = *this,
                    .anchor = anchor,
                    .pex_node = pex_node,
                    .query = query,
                    .orientation = orientation,
                    .reference = references[reference_id],
                    .already_verified_intervals = already_verified_intervals_per_reference[reference_id],
                    .alignments = alignments,
                    .stats = stats
                };

                verifier.verify(config.verification_kind);
            }
        }
    }
}

} // namespace pex
