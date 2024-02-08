#include <pex.hpp>

#include <iostream>
#include <utility>
#include <sstream>

size_t ceil_div(size_t const a, size_t const b) {
    return (a % b) ? a / b + 1 : a / b;
}

pex_tree::pex_tree(pex_tree_config const config) 
    : leaf_query_length{config.total_query_length / (config.query_num_errors + 1)},
    leaf_num_errors{config.leaf_num_errors} {
    // use 1 based indices until final computation to make sure to match pseudocode
    add_nodes(
        1, 
        config.total_query_length,
        config.query_num_errors,
        null_id
    );
}

std::string pex_tree::node::to_string() const {
    std::stringstream s{};
    s << "{ parent_id: " << parent_id <<
        ", from: " << query_index_from <<
        ", to: " << query_index_to <<
        ", errors: " << num_errors << " }";
    return s.str();
}


void pex_tree::debug_print() const{
    std::cout << "--- INNER NODES: ---\n";
    for (auto const& node : inner_nodes) {
        std::cout << node.to_string() << '\n';
    }

    std::cout << "--- LEAF NODES: ---\n";
    for (auto const& node : leaves) {
        std::cout << node.to_string() << '\n';
    }
}

std::vector<pex_tree::node> const& pex_tree::get_leaves() const {
    return leaves;
}

void pex_tree::add_nodes(
    size_t const query_index_from,
    size_t const query_index_to,
    size_t const num_errors, 
    size_t const parent_id
) {
    // is this the correct meaning of this variable?
    size_t const num_leafs_left = ceil_div(num_errors + 1, 2);

    node const curr_node = {
        parent_id,
        query_index_from - 1, // transform to 0-based index
        query_index_to - 1, // transform to 0-based index
        num_errors
    };

    if (num_errors <= leaf_num_errors) {
        leaves.push_back(curr_node);
    } else {
        size_t const curr_node_id = inner_nodes.size();
        inner_nodes.push_back(curr_node);

        size_t const query_split_index = query_index_from + num_leafs_left * leaf_query_length;

        add_nodes(
            query_index_from,
            query_split_index - 1,
            (num_leafs_left * num_errors) / (num_errors + 1),
            curr_node_id
        );
        add_nodes(
            query_split_index,
            query_index_to,
            ((num_errors + 1 - num_leafs_left) * num_errors) / (num_errors + 1),
            curr_node_id
        );
    }
}

pex_tree const& pex_tree_cache::get(pex_tree_config const config) {
    auto [iter, _] = trees.try_emplace(config.total_query_length, config);
    return iter->second;
}
