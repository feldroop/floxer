#include <alignment_algorithm.hpp>
#include <pex.hpp>

#include <cassert>
#include <ranges>
#include <tuple>

#include <spdlog/spdlog.h>

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

bool pex_tree::node::is_root() const {
    return parent_id == null_id;
}

void pex_tree::search(
    std::vector<input::reference_record> const& references,
    std::span<const uint8_t> const fastq_query,
    alignment::fastq_query_alignments& output_alignments,
    bool const is_reverse_complement,
    search::search_scheme_cache& scheme_cache,
    fmindex const& index
) const {
    auto const leaf_queries = generate_leaf_queries(fastq_query);

    spdlog::trace("searching seeds in FM-index");

    auto const hits = search::search_leaf_queries(
        leaf_queries,
        index,
        scheme_cache,
        references.size()
    );

    size_t num_hits = 0;

    for (auto const& leaf_query_hits : hits) {
        for (auto const& leaf_to_reference_hits : leaf_query_hits) {
            num_hits += leaf_to_reference_hits.size();
        }
    }

    spdlog::trace("found {} hits, now verifying/aligning", num_hits);

    for (size_t leaf_query_id = 0; leaf_query_id < leaf_queries.size(); ++leaf_query_id) {
        for (size_t reference_id = 0; reference_id < references.size(); ++reference_id) {
            for (auto const& hit : hits[leaf_query_id][reference_id]) {
                hierarchical_verification(
                    hit,
                    leaf_query_id,
                    fastq_query,
                    references[reference_id],
                    output_alignments,
                    is_reverse_complement
                );
            }
        }
    }
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

std::tuple<size_t, size_t> compute_reference_span_start_and_length(
    search::hit const& hit,
    pex_tree::node const& pex_node,
    size_t const leaf_query_index_from,
    size_t const full_reference_length
) {
    // this extra wiggle room is added around the reference span because it was observed
    // that it leads no nicer alignments in certain edge cases 
    // and it likely has no impact on performance
    size_t constexpr extra_wiggle_room = 1;

    int64_t const start_signed = static_cast<int64_t>(hit.position) - 
        (leaf_query_index_from - pex_node.query_index_from)
        - pex_node.num_errors - extra_wiggle_room;
    size_t const reference_span_start = start_signed >= 0 ? start_signed : 0;
    size_t const reference_span_length = std::min(
        pex_node.query_length() + 2 * pex_node.num_errors + 1 + 2 * extra_wiggle_room,
        full_reference_length - reference_span_start
    );

    return std::make_tuple(reference_span_start, reference_span_length);
}

// returns whether an alignment was found
bool try_to_align_corresponding_query_span_at_anchor(
    search::hit const& hit,
    pex_tree::node const& pex_node,
    size_t const leaf_query_index_from,
    input::reference_record const& reference,
    std::span<const uint8_t> const fastq_query,
    alignment::fastq_query_alignments& alignments,
    bool const is_reverse_complement
) {
    auto const full_reference_span = std::span<const uint8_t>(reference.rank_sequence);

    auto const [reference_span_start, reference_span_length] = compute_reference_span_start_and_length(
        hit,
        pex_node,
        leaf_query_index_from,
        full_reference_span.size()
    );

    auto const this_node_query_span = fastq_query.subspan(
        pex_node.query_index_from,
        pex_node.query_length()
    );

    auto const reference_subspan =
        full_reference_span.subspan(reference_span_start, reference_span_length);

    auto alignment_insertion_gatekeeper = alignments.get_insertion_gatekeeper(
        reference.id,
        reference_span_start,
        reference_span_length,
        is_reverse_complement
    );

    bool const query_found = alignment::VerifyingAligner::align_query(
        // the reversing here is done to allow the DP traceback to start from the start position
        // of the alignment in the reference. This in turn allows to skip tracebacks for many
        // alignment candidates that are not useful
        std::views::reverse(reference_subspan),
        std::views::reverse(this_node_query_span),
        pex_node.num_errors,
        pex_node.is_root(),
        alignment_insertion_gatekeeper
    );

    return query_found;
}

void pex_tree::hierarchical_verification(
    search::hit const& hit,
    size_t const leaf_query_id,
    std::span<const uint8_t> const fastq_query,
    input::reference_record const& reference,
    alignment::fastq_query_alignments& alignments,
    bool const is_reverse_complement
) const {
    // this depends on the implementation of generate_leave_queries returning the
    // leaf queries in the same order as the leaves (which it should always do!)
    auto pex_node = leaves.at(leaf_query_id);
    size_t const leaf_query_index_from = pex_node.query_index_from;

    // case for when the whole PEX tree is just a single root
    // this could be optimized by not aligning again, but instead using the FM-index alignment
    if (pex_node.is_root()) {
        bool const query_found = try_to_align_corresponding_query_span_at_anchor(
            hit,
            pex_node,
            leaf_query_index_from,
            reference,
            fastq_query,
            alignments,
            is_reverse_complement
        );

        (void)query_found;
        assert(query_found);

        return;
    }

    pex_node = inner_nodes.at(pex_node.parent_id);

    while (true) {
        bool const query_found = try_to_align_corresponding_query_span_at_anchor(
            hit,
            pex_node,
            leaf_query_index_from,
            reference,
            fastq_query,
            alignments,
            is_reverse_complement
        );

        if (!query_found || pex_node.is_root()) {
            break;
        }

        pex_node = inner_nodes.at(pex_node.parent_id);
    }
}

pex_tree const& pex_tree_cache::get(pex_tree_config const config) {
    auto [iter, _] = trees.try_emplace(config.total_query_length, config);
    return iter->second;
}
