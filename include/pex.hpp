#pragma once

#include <fmindex.hpp>
#include <io.hpp>
#include <search.hpp>
#include <verification.hpp>

#include <limits>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

#include <fmindex-collection/concepts.h>

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
        size_t parent_id;
        size_t query_index_from;
        size_t query_index_to;
        size_t num_errors;

        size_t query_length() const;
    };

    // return[reference_id][alignment_end_position] -> alignment data
    std::vector<std::unordered_map<size_t, verification::query_alignment>> search(
        std::vector<io::record> const& references,
        std::span<const uint8_t> const fastq_query,
        search::search_scheme_cache& scheme_cache,
        fmindex const& index
    ) const;

private:
    static constexpr size_t null_id = std::numeric_limits<size_t>::max();
    
    std::vector<node> inner_nodes;
    std::vector<node> leaves;

    // this refers to the original version where leaves have 0 errors 
    size_t const no_error_leaf_query_length;
    // this is the leaf query length we get with the version where
    // the leaves might have > 0 errors
    size_t actual_leaf_query_length;
    
    size_t const leaf_num_errors;

    void add_nodes(
        size_t const query_index_from,
        size_t const query_index_to,
        size_t const num_errors, 
        size_t const parent_id
    );

    // returns queries in the same order as the leavesare stored in the tree
    std::vector<std::span<const uint8_t>> generate_leaf_queries(std::span<const uint8_t> const& full_query) const;

    void hierarchical_verification(
        search::hit const& hit,
        size_t const leaf_query_id,
        std::span<const uint8_t> const fastq_query,
        std::span<const uint8_t> const reference,
        std::unordered_map<size_t, verification::query_alignment>& reference_alignments
    ) const;
};

class pex_tree_cache {
public:
    pex_tree const& get(pex_tree_config const config);
private:
    // query length determines PEX tree structure uniquely in this application
    std::unordered_map<size_t, pex_tree> trees;
};
