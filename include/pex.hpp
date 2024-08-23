#pragma once

#include <alignment.hpp>
#include <fmindex.hpp>
#include <input.hpp>
#include <intervals.hpp>
#include <search.hpp>
#include <statistics.hpp>

#include <limits>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

#include <fmindex-collection/concepts.h>

namespace pex {

enum pex_tree_build_strategy {
    recursive,
    bottom_up
};

struct pex_tree_config {
    size_t const total_query_length;
    size_t const query_num_errors;
    size_t const leaf_max_num_errors;
    pex_tree_build_strategy const build_strategy;
};

enum class verification_kind_t {
    direct_full, hierarchical
};

struct pex_alignment_config {
    search::searcher const& searcher;
    intervals::use_interval_optimization const use_interval_optimization;
    verification_kind_t const verification_kind;
    double const extra_verification_ratio;
};

// based on chapter 6.5.1 from the book "Flexible Pattern Matching in Strings" by Navarro and Raffinot
// DOI: https://doi.org/10.1017/CBO9781316135228
class pex_tree {
public:
    struct node {
        static constexpr size_t null_id = std::numeric_limits<size_t>::max();

        // parent_id is index of parent node in inner_nodes or null_id for the root
        size_t parent_id;

        // incluvsive range [from, to]
        size_t query_index_from;
        size_t query_index_to;

        // for search with FM-index
        size_t num_errors;

        size_t length_of_query_span() const;
        bool is_root() const;

        std::string dot_statement(size_t const id) const;
    };

    pex_tree(pex_tree_config const config);

    node const& root() const;

    // throws std::runtime_error if root is given
    node const& get_parent_of_child(node const& child) const;

    std::vector<node> const& get_leaves() const;

    // returns seeds in the same order as the leaves are stored in the tree (index in vector = seed_id)
    std::vector<search::seed> generate_seeds(std::span<const uint8_t> const query) const;

    alignment::query_alignments align_forward_and_reverse_complement(
        std::vector<input::reference_record> const& references,
        std::span<const uint8_t> const query,
        pex_alignment_config const config,
        statistics::search_and_alignment_statistics& stats
    ) const;

    std::string dot_statement() const;
private:

    size_t num_leaves() const;

    void add_node_to_dot_statement(node const& curr_node, size_t const id, std::string& dot) const;

    // this is the building from the book with a small adjustment for varying number of leaf errors
    void add_nodes_recursive(
        size_t const query_index_from,
        size_t const query_index_to,
        size_t const num_errors,
        size_t const parent_id
    );

    void add_nodes_bottom_up(pex_tree_config const& config);

    // for bottom_up build strategy
    void create_leaves(
        pex_tree_config const& config,
        size_t const num_desired_leaves,
        size_t const num_leaves_with_less_errors
    );

    // for bottom up build strategy, returns parent node for nodes in child_nodes and sets their parent_id
    node create_parent_node(std::span<node> const child_nodes, size_t const parent_id);

    void align_query_in_given_orientation(
        std::vector<input::reference_record> const& references,
        std::span<const uint8_t> const query,
        alignment::query_alignments& alignments,
        alignment::query_orientation const orientation,
        pex_alignment_config const config,
        statistics::search_and_alignment_statistics& stats
    ) const;

    std::vector<node> inner_nodes;
    std::vector<node> leaves;

    // this refers to the original version where leaves have 0 errors
    size_t const no_error_seed_length;

    size_t const leaf_max_num_errors;
};

class pex_tree_cache {
public:
    pex_tree const& get(pex_tree_config const config);
private:
    // query length determines PEX tree structure uniquely in this application
    // because either there is a constant number of errors, or the number
    // of errors per query is a function of only the query length
    std::unordered_map<size_t, pex_tree> trees;
};

} // namespace pex
