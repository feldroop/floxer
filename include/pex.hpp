#pragma once

#include <alignment.hpp>
#include <fmindex.hpp>
#include <input.hpp>
#include <search.hpp>
#include <statistics.hpp>

#include <limits>
#include <map>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

#include <fmindex-collection/concepts.h>

namespace pex {

struct pex_tree_config {
    size_t const total_query_length;
    size_t const query_num_errors;
    size_t const leaf_max_num_errors;
};

class pex_tree {
public:
    pex_tree() = delete;
    pex_tree(pex_tree_config const config);

    struct node {
        size_t parent_id;
        size_t query_index_from;
        size_t query_index_to;
        size_t num_errors;

        size_t length_of_query_span() const;
        bool is_root() const;

        std::string dot_statement(size_t const id) const;
    };

    alignment::query_alignments align_forward_and_reverse_complement(
        std::vector<input::reference_record> const& references,
        std::span<const uint8_t> const query,
        search::searcher const& searcher,
        statistics::search_and_alignment_statistics& stats
    ) const;

    std::string dot_statement() const;

    size_t max_leaf_query_span() const;

    size_t num_leaves() const;

private:
    static constexpr size_t null_id = std::numeric_limits<size_t>::max();
    
    std::vector<node> inner_nodes;
    std::vector<node> leaves;

    // this refers to the original version where leaves have 0 errors 
    size_t const no_error_leaf_query_length;
    
    size_t const leaf_max_num_errors;

    void add_nodes(
        size_t const query_index_from,
        size_t const query_index_to,
        size_t const num_errors, 
        size_t const parent_id
    );

    void add_node_to_dot_statement(node const& curr_node, size_t const id, std::string& dot) const;

    // returns seeds in the same order as the leaves are stored in the tree
    std::vector<search::seed> generate_seeds(std::span<const uint8_t> const query) const;

    void align_query_in_given_orientation(
        std::vector<input::reference_record> const& references,
        std::span<const uint8_t> const query,
        alignment::query_alignments& alignments,
        bool const is_reverse_complement,
        search::searcher const& searcher,
        statistics::search_and_alignment_statistics& stats
    ) const;

    void hierarchical_verification(
        search::anchor const& anchor,
        size_t const seed_id,
        std::span<const uint8_t> const query,
        input::reference_record const& reference,
        alignment::query_alignments& alignments,
        bool const is_reverse_complement
    ) const;
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

size_t ceil_div(size_t const a, size_t const b);

struct span_config {
    size_t const offset{};
    size_t const length{};
};

span_config compute_reference_span_start_and_length(
    search::anchor const& anchor,
    pex_tree::node const& pex_node,
    size_t const leaf_query_index_from,
    size_t const full_reference_length
);

// returns whether an alignment was found
bool try_to_align_corresponding_query_span_at_anchor(
    search::anchor const& anchor,
    pex_tree::node const& pex_node,
    size_t const seed_query_index_from,
    input::reference_record const& reference,
    std::span<const uint8_t> const query,
    alignment::query_alignments& alignments,
    bool const is_reverse_complement
);

} // namespace internal

} // namespace pex
