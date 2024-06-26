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

// based on chapter 6.5.1 from the book "Flexible Pattern Matching in Strings" by Navarro and Raffinot
// DOI: https://doi.org/10.1017/CBO9781316135228
class pex_tree {
public:
    struct node {
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

    alignment::query_alignments align_forward_and_reverse_complement(
        std::vector<input::reference_record> const& references,
        std::span<const uint8_t> const query,
        search::searcher const& searcher,
        intervals::use_interval_optimization const use_interval_optimization,
        statistics::search_and_alignment_statistics& stats
    ) const;

    std::string dot_statement() const;

private:
    node const& root() const;

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

    // returns seeds in the same order as the leaves are stored in the tree (index in vector = seed_id)
    std::vector<search::seed> generate_seeds(std::span<const uint8_t> const query) const;

    void align_query_in_given_orientation(
        std::vector<input::reference_record> const& references,
        std::span<const uint8_t> const query,
        alignment::query_alignments& alignments,
        alignment::query_orientation const orientation,
        search::searcher const& searcher,
        intervals::use_interval_optimization const use_interval_optimization,
        statistics::search_and_alignment_statistics& stats
    ) const;

    void hierarchical_verification(
        search::anchor_t const& anchor,
        size_t const seed_id,
        std::span<const uint8_t> const query,
        alignment::query_orientation const orientation,
        input::reference_record const& reference,
        intervals::verified_intervals& already_verified_intervals,
        alignment::query_alignments& alignments,
        statistics::search_and_alignment_statistics& stats
    ) const;

    static constexpr size_t null_id = std::numeric_limits<size_t>::max();

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

namespace internal {

struct span_config {
    size_t const offset{};
    size_t const length{};

    intervals::half_open_interval as_half_open_interval() const;
};

span_config compute_reference_span_start_and_length(
    search::anchor_t const& anchor,
    pex_tree::node const& pex_node,
    size_t const leaf_query_index_from,
    size_t const full_reference_length,
    size_t const extra_wiggle_room
);

alignment::alignment_outcome try_to_align_pex_node_query_with_reference_span(
    pex_tree::node const& pex_node,
    input::reference_record const& reference,
    span_config const reference_span_config,
    std::span<const uint8_t> const query,
    alignment::query_orientation const orientation,
    alignment::query_alignments& alignments,
    statistics::search_and_alignment_statistics& stats
);

} // namespace internal

} // namespace pex
