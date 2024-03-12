#include <pex.hpp>

#include <fmt/core.h>

size_t ceil_div(size_t const a, size_t const b) {
    return (a % b) ? a / b + 1 : a / b;
}

pex_tree::pex_tree(pex_tree_config const config) 
    : no_error_leaf_query_length{config.total_query_length / (config.query_num_errors + 1)},
    leaf_max_num_errors{config.leaf_max_num_errors} {
    // use 1 based indices until final computation to make sure to match pseudocode
    add_nodes(
        1, 
        config.total_query_length,
        config.query_num_errors,
        null_id
    );
}

size_t pex_tree::node::query_length() const {
    return query_index_to - query_index_from + 1;
}

std::vector<std::map<size_t, verification::query_alignment>> pex_tree::search(
    std::vector<io::record> const& references,
    std::span<const uint8_t> const fastq_query,
    search::search_scheme_cache& scheme_cache,
    fmindex const& index
) const {
    auto const leaf_queries = generate_leaf_queries(fastq_query);
    
    auto const hits = search::search_leaf_queries(
        leaf_queries,
        index,
        scheme_cache,
        references.size()
    );

    // alignments[reference_id][end_position] -> alignment of fastq query to this reference
    std::vector<std::map<size_t, verification::query_alignment>> alignments(
        references.size()
    );

    for (size_t leaf_query_id = 0; leaf_query_id < leaf_queries.size(); ++leaf_query_id) {
        for (size_t reference_id = 0; reference_id < references.size(); ++reference_id) {
            auto const reference = std::span<const uint8_t>(references[reference_id].sequence);
            auto & references_alignments = alignments[reference_id];

            for (auto const& hit : hits[leaf_query_id][reference_id]) {
                hierarchical_verification(
                    hit,
                    leaf_query_id,
                    fastq_query,
                    reference,
                    references_alignments
                );
            }
        }
    }

    return alignments;
}

void pex_tree::add_nodes(
    size_t const query_index_from,
    size_t const query_index_to,
    size_t const num_errors, 
    size_t const parent_id
) {
    // not sure that this name is the correct meaning of this value from the book
    size_t const num_leafs_left = ceil_div(num_errors + 1, 2);

    node const curr_node = {
        parent_id,
        query_index_from - 1, // transform to 0-based index
        query_index_to - 1, // transform to 0-based index
        num_errors
    };

    if (num_errors <= leaf_max_num_errors) {
        leaves.push_back(curr_node);
    } else {
        size_t const curr_node_id = inner_nodes.size();
        inner_nodes.push_back(curr_node);

        size_t const query_split_index = query_index_from + num_leafs_left * no_error_leaf_query_length;

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

std::vector<search::query> pex_tree::generate_leaf_queries(
    std::span<const uint8_t> const& full_query
) const {
    std::vector<search::query> leaf_queries{};
    leaf_queries.reserve(leaves.size());

    for (auto const& leaf : leaves) {
        auto const leaf_query_span = full_query.subspan(leaf.query_index_from, leaf.query_length());
        leaf_queries.emplace_back(std::move(leaf_query_span), leaf.num_errors);
    }   

    return leaf_queries;
}

void pex_tree::hierarchical_verification(
    search::hit const& hit,
    size_t const leaf_query_id,
    std::span<const uint8_t> const fastq_query,
    std::span<const uint8_t> const reference,
    std::map<size_t, verification::query_alignment>& reference_alignments
) const {    
    // this depends on the implementation of generate_leave_queries returning the
    // leaf queries in the same order as the leaves (which it should always do!)
    auto pex_node = leaves.at(leaf_query_id);
    size_t const leaf_query_index_from = pex_node.query_index_from;
    pex_node = inner_nodes[pex_node.parent_id];

    while (true) {
        int64_t const start_signed = static_cast<int64_t>(hit.position) - 
            (leaf_query_index_from - pex_node.query_index_from)
            - pex_node.num_errors;
        size_t const reference_span_start = start_signed >= 0 ? start_signed : 0;
        size_t const reference_span_length = std::min(
            pex_node.query_length() + 2 * pex_node.num_errors + 1,
            reference.size() - reference_span_start
        );
        auto const& this_node_query = fastq_query.subspan(
            pex_node.query_index_from,
            pex_node.query_length()
        );

        bool const curr_node_is_root = pex_node.parent_id == null_id;

        auto alignments_wrapper = verification::alignment_output_gatekeeper(
            reference_span_start, reference_alignments
        );

        bool const query_found = verification::align_query(
            reference.subspan(reference_span_start, reference_span_length),
            this_node_query,
            pex_node.num_errors,
            curr_node_is_root,
            alignments_wrapper // useful alignments are written into this
        );

        if (!query_found || curr_node_is_root) {
            break;
        }

        pex_node = inner_nodes.at(pex_node.parent_id);
    }
}

pex_tree const& pex_tree_cache::get(pex_tree_config const config) {
    auto [iter, _] = trees.try_emplace(config.total_query_length, config);
    return iter->second;
}
