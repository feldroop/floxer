#pragma once

#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

struct pex_tree_config {
    size_t const total_query_length;
    size_t const query_num_errors;
    size_t const leaf_num_errors = 0;
};

class pex_tree {
public:
    pex_tree() = delete;
    pex_tree(pex_tree_config const config);

    struct node {
        size_t const parent_id;
        size_t const query_index_from;
        size_t const query_index_to;
        size_t const num_errors;

        std::string to_string() const;
    };

    void debug_print() const;

    std::vector<node> const& get_leaves() const;

private:
    static constexpr size_t null_id = std::numeric_limits<size_t>::max();
    
    std::vector<node> inner_nodes;
    std::vector<node> leaves;

    size_t const leaf_query_length;
    size_t const leaf_num_errors;

    void add_nodes(
        size_t const query_index_from,
        size_t const query_index_to,
        size_t const num_errors, 
        size_t const parent_id
    );
};

class pex_tree_cache {
public:
    pex_tree const& get(pex_tree_config const config);
private:
    // query length determines PEX tree structure uniquely in this application
    std::unordered_map<size_t, pex_tree> trees;
};
