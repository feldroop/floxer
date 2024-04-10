#pragma once

#include <fmindex.hpp>
#include <input.hpp>
#include <search.hpp>
#include <alignment.hpp>

#include <limits>
#include <map>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

#include <fmindex-collection/concepts.h>

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

        size_t query_length() const;
        bool is_root() const;
    };

    void search(
        std::vector<input::reference_record> const& references,
        std::span<const uint8_t> const fastq_query,
        alignment::fastq_query_alignments& output_alignments,
        bool const is_reverse_complement,
        search::search_scheme_cache& scheme_cache,
        fmindex const& index
    ) const;

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

    // returns queries in the same order as the leaves are stored in the tree
    std::vector<search::query> generate_leaf_queries(std::span<const uint8_t> const& full_query) const;

    void hierarchical_verification(
        search::hit const& hit,
        size_t const leaf_query_id,
        std::span<const uint8_t> const fastq_query,
        input::reference_record const& reference,
        alignment::fastq_query_alignments& alignments,
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
