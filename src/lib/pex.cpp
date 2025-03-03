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

pex_tree_config::pex_tree_config(size_t const query_sequence_length, cli::command_line_input const& cli_input)
    : total_query_length{query_sequence_length},
        query_num_errors{input::num_errors_from_user_config(query_sequence_length, cli_input)},
        leaf_max_num_errors{cli_input.pex_seed_num_errors()},
        build_strategy{
            cli_input.bottom_up_pex_tree_building() ?
                pex::pex_tree_build_strategy::bottom_up :
                pex::pex_tree_build_strategy::recursive
        }
{}

pex_tree_config::pex_tree_config(
    size_t const total_query_length_,
    size_t const query_num_errors_,
    size_t const leaf_max_num_errors_,
    pex_tree_build_strategy const build_strategy_
) : total_query_length{total_query_length_},
    query_num_errors{query_num_errors_},
    leaf_max_num_errors{leaf_max_num_errors_},
    build_strategy{build_strategy_}
{}


pex_verification_config::pex_verification_config(cli::command_line_input const& cli_input)
    : use_interval_optimization{
        cli_input.use_interval_optimization() ?
            intervals::use_interval_optimization::on :
            intervals::use_interval_optimization::off
    },
    verification_kind{
        cli_input.direct_full_verification() ?
            pex::verification_kind_t::direct_full :
            pex::verification_kind_t::hierarchical
    },
    extra_verification_ratio{cli_input.extra_verification_ratio()}
{}

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

std::vector<pex_tree::node> const& pex_tree::get_leaves() const {
    return leaves;
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
                node::null_id
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
        leaves.emplace_back(node {
            .parent_id = node::null_id,
            .query_index_from = 0,
            .query_index_to = config.total_query_length - 1,
            .num_errors = config.query_num_errors
        });

        return;
    }

    create_leaves(config, num_desired_leaves);

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
    inner_nodes.front().parent_id = node::null_id;

    // free reserved space that is not needed
    inner_nodes.shrink_to_fit();
}

void pex_tree::create_leaves(
    pex_tree_config const& config,
    size_t const num_desired_leaves
) {
    size_t const base_seed_length = config.total_query_length / num_desired_leaves;
    size_t const seed_length_remainder = config.total_query_length % num_desired_leaves;

    leaves.reserve(num_desired_leaves);

    size_t current_start_index = 0;
    for (size_t i = 0; i < num_desired_leaves; ++i) {
        size_t const curr_leaf_length = i < seed_length_remainder ?
            base_seed_length + 1 :
            base_seed_length;

        leaves.emplace_back(node{
            .parent_id = 0, // will be set later
            .query_index_from = current_start_index,
            .query_index_to = current_start_index + curr_leaf_length - 1,
            .num_errors = config.leaf_max_num_errors
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
    std::span<const uint8_t> const query,
    size_t const seed_sampling_step_size
) const {
    std::vector<search::seed> seeds{};
    seeds.reserve(leaves.size());

    for (size_t pex_leaf_index = 0; pex_leaf_index < leaves.size(); pex_leaf_index += seed_sampling_step_size) {
        auto const& leaf = leaves[pex_leaf_index];
        auto const seed_span = query.subspan(leaf.query_index_from, leaf.length_of_query_span());
        seeds.emplace_back(search::seed {
            .sequence = seed_span,
            .num_errors = leaf.num_errors,
            .query_position = leaf.query_index_from,
            .pex_leaf_index = pex_leaf_index
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

} // namespace pex
