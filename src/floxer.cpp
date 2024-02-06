#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <vector>

// #include <fmindex-collection.h>

size_t ceil_div(size_t const a, size_t const b) {
    return (a % b) ? a / b + 1 : a / b;
}

class pex_tree {
public:
    pex_tree() = delete;
    pex_tree(
        size_t const total_query_length,
        uint8_t const num_query_errors
        // uint8_t const num_errors_for_search = 0
    ) : leaf_query_length{total_query_length / (num_query_errors + 1)} {
        // use 1 based indices until final computation to make sure to match pseudocode
        add_nodes(
            1, 
            total_query_length,
            num_query_errors,
            null_id
        );
    }

    struct node {
        size_t const parent_id;
        size_t const query_index_from;
        size_t const query_index_to;
        uint8_t const num_errors;

        std::string to_string() const {
            std::stringstream s{};
            s << "{ parent_id: " << parent_id <<
                ", from: " << query_index_from <<
                ", to: " << query_index_to <<
                ", errors: " << static_cast<int>(num_errors) << " }";
            return s.str();
        }
    };

    void debug_print() {
        std::cout << "--- INNER NODES: ---\n";
        for (auto const& node : inner_nodes) {
            std::cout << node.to_string() << '\n';
        }

        std::cout << "--- LEAF NODES: ---\n";
        for (auto const& node : leafs) {
            std::cout << node.to_string() << '\n';
        }
    }

private:
    static constexpr size_t null_id = std::numeric_limits<size_t>::max();
    
    std::vector<node> inner_nodes;
    std::vector<node> leafs;

    size_t const leaf_query_length; 

    void add_nodes(
        size_t const query_index_from,
        size_t const query_index_to,
        uint8_t const num_errors, 
        size_t const parent_id
    ) {
        // is this the correct meaning of this variable?
        size_t const num_leafs_left = ceil_div(static_cast<size_t>(num_errors + 1), 2);

        node const curr_node = {
            parent_id,
            query_index_from - 1, // transform to 0-based index
            query_index_to - 1, // transform to 0-based index
            num_errors
        };

        if (num_errors == 0) {
            leafs.push_back(curr_node);
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
};
