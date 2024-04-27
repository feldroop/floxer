#include <alignment_algorithm.hpp>
#include <pex.hpp>

#include <cassert>
#include <ranges>

#include <ivsigma/ivsigma.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/spdlog.h>

namespace pex {

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

size_t pex_tree::node::length_of_query_span() const {
    return query_index_to - query_index_from + 1;
}

bool pex_tree::node::is_root() const {
    return parent_id == null_id;
}

std::string pex_tree::node::dot_statement(size_t const id) const {
    return fmt::format(
        "{} [label=\"errors: {}\\nlength: {}\\nrange: [{},{}]\"];\n",
        id,
        num_errors,
        length_of_query_span(),
        query_index_from,
        query_index_to
    );
}

alignment::query_alignments pex_tree::align_forward_and_reverse_complement(
    std::vector<input::reference_record> const& references,
    std::span<const uint8_t> const query,
    search::searcher const& searcher,
    statistics::search_and_alignment_statistics& stats
) const {
    auto alignments = alignment::query_alignments(references.size());

    bool is_reverse_complement = false;
    align_query_in_given_orientation(
        references,
        query,
        alignments,
        is_reverse_complement,
        searcher,
        stats
    );

    auto const reverse_complement_query = 
        ivs::reverse_complement_rank<ivs::d_dna4>(query);
    is_reverse_complement = true;

    align_query_in_given_orientation(
        references,
        reverse_complement_query,
        alignments,
        is_reverse_complement,
        searcher,
        stats
    );

    stats.add_num_alignments(alignments.size());

    for (size_t reference_id = 0; reference_id < references.size(); ++reference_id) {
        for (auto const& [_, alignment] : alignments.to_reference(reference_id)) {
            stats.add_alignment_edit_distance(alignment.num_errors);
        }
    }

    return alignments;
}

std::string pex_tree::dot_statement() const {
    std::string dot = fmt::format(
        "graph {{\n"
        "label = \"PEX tree for query length {}, {} errors and leaf threshold {} ({} leaves)\";\n"
        "labelloc = \"t\";\n"
        "node [shape=record];\n",
        inner_nodes[0].query_index_to + 1,
        inner_nodes[0].num_errors,
        leaf_max_num_errors,
        num_leaves()
    );
    
    size_t id = 0;
    for (auto const& inner_node : inner_nodes) {
        add_node_to_dot_statement(inner_node, id, dot);
        ++id;
    }
    for (auto const& leaf_node : leaves) {
        add_node_to_dot_statement(leaf_node, id, dot);
        ++id;
    }

    dot += "}\n";

    return dot;
}

void pex_tree::align_query_in_given_orientation(
    std::vector<input::reference_record> const& references,
    std::span<const uint8_t> const query,
    alignment::query_alignments& alignments,
    bool const is_reverse_complement,
    search::searcher const& searcher,
    statistics::search_and_alignment_statistics& stats
) const {
    auto const seeds = generate_seeds(query);
    stats.add_seed_lengths_and_num_seeds_per_query(seeds);

    auto const search_result = searcher.search_seeds(seeds);
    stats.add_statistics_for_search_result(search_result);

    for (size_t seed_id = 0; seed_id < seeds.size(); ++seed_id) {
        auto const& anchors_of_seed = search_result.anchors_by_seed[seed_id];
        
        if (anchors_of_seed.excluded) {
            continue;
        }
        
        for (size_t reference_id = 0; reference_id < references.size(); ++reference_id) {
            for (auto const& anchor : anchors_of_seed.anchors_by_reference[reference_id]) {
                hierarchical_verification(
                    anchor,
                    seed_id,
                    query,
                    references[reference_id],
                    alignments,
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
    size_t const num_leafs_left = internal::ceil_div(num_errors + 1, 2);

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

void pex_tree::add_node_to_dot_statement(node const& curr_node, size_t const id, std::string& dot) const {
    dot += curr_node.dot_statement(id);
    if (!curr_node.is_root()) {
        dot += fmt::format("{} -- {};\n", id, curr_node.parent_id);
    }
}

std::vector<search::seed> pex_tree::generate_seeds(
    std::span<const uint8_t> const query
) const {
    std::vector<search::seed> seeds{};
    seeds.reserve(leaves.size());

    for (auto const& leaf : leaves) {
        auto const seed_span = query.subspan(leaf.query_index_from, leaf.length_of_query_span());
        seeds.emplace_back(seed_span, leaf.num_errors);
    }

    return seeds;
}

void pex_tree::hierarchical_verification(
    search::anchor const& anchor,
    size_t const seed_id,
    std::span<const uint8_t> const query,
    input::reference_record const& reference,
    alignment::query_alignments& alignments,
    bool const is_reverse_complement
) const {
    // this depends on the implementation of generate_leave_queries returning the
    // leaf queries in the same order as the leaves (which it should always do!)
    auto pex_node = leaves.at(seed_id);
    size_t const seed_query_index_from = pex_node.query_index_from;

    // case for when the whole PEX tree is just a single root
    // this could be optimized by not aligning again, but instead using the FM-index alignment
    if (pex_node.is_root()) {
        [[maybe_unused]] bool const query_found = internal::try_to_align_corresponding_query_span_at_anchor(
            anchor,
            pex_node,
            seed_query_index_from,
            reference,
            query,
            alignments,
            is_reverse_complement
        );

        assert(query_found);

        return;
    }

    pex_node = inner_nodes.at(pex_node.parent_id);

    while (true) {
        bool const query_found = internal::try_to_align_corresponding_query_span_at_anchor(
            anchor,
            pex_node,
            seed_query_index_from,
            reference,
            query,
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

namespace internal {

size_t ceil_div(size_t const a, size_t const b) {
    return (a % b) ? a / b + 1 : a / b;
}

span_config compute_reference_span_start_and_length(
    search::anchor const& anchor,
    pex_tree::node const& pex_node,
    size_t const leaf_query_index_from,
    size_t const full_reference_length
) {
    // this extra wiggle room is added around the reference span because it was observed
    // that it leads no nicer alignments in certain edge cases 
    // and it likely has no impact on performance
    size_t constexpr extra_wiggle_room = 1;

    int64_t const start_signed = static_cast<int64_t>(anchor.position) - 
        (leaf_query_index_from - pex_node.query_index_from)
        - pex_node.num_errors - extra_wiggle_room;
    size_t const reference_span_start = start_signed >= 0 ? start_signed : 0;
    size_t const reference_span_length = std::min(
        pex_node.length_of_query_span() + 2 * pex_node.num_errors + 1 + 2 * extra_wiggle_room,
        full_reference_length - reference_span_start
    );

    return span_config{
        .offset = reference_span_start,
        .length = reference_span_length
    };
}

bool try_to_align_corresponding_query_span_at_anchor(
    search::anchor const& anchor,
    pex_tree::node const& pex_node,
    size_t const seed_query_index_from,
    input::reference_record const& reference,
    std::span<const uint8_t> const query,
    alignment::query_alignments& alignments,
    bool const is_reverse_complement
) {
    auto const full_reference_span = std::span<const uint8_t>(reference.rank_sequence);

    auto const reference_span_config = compute_reference_span_start_and_length(
        anchor,
        pex_node,
        seed_query_index_from,
        full_reference_span.size()
    );

    auto const this_node_query_span = query.subspan(
        pex_node.query_index_from,
        pex_node.length_of_query_span()
    );

    auto const reference_subspan = full_reference_span.subspan(
        reference_span_config.offset,
        reference_span_config.length
    );

    auto alignment_insertion_gatekeeper = alignments.get_insertion_gatekeeper(
        reference.id,
        reference_span_config.offset,
        reference_span_config.length,
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

} // namespace internal

} // namespace pex
