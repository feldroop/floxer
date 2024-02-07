#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <vector>

class pex_tree {
public:
    pex_tree() = delete;
    pex_tree(
        size_t const total_query_length,
        size_t const num_query_errors,
        size_t const leaf_num_errors_ = 0
    );

    struct node {
        size_t const parent_id;
        size_t const query_index_from;
        size_t const query_index_to;
        size_t const num_errors;

        std::string to_string() const;
    };

    void debug_print();

private:
    static constexpr size_t null_id = std::numeric_limits<size_t>::max();
    
    std::vector<node> inner_nodes;
    std::vector<node> leafs;

    size_t const leaf_query_length;
    size_t const leaf_num_errors;

    void add_nodes(
        size_t const query_index_from,
        size_t const query_index_to,
        size_t const num_errors, 
        size_t const parent_id
    );
};
