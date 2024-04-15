#include <fmindex.hpp>
#include <pex.hpp>
#include <search.hpp>

#include <algorithm>
#include <limits>
#include <ranges>

#include <fmindex-collection/search/SearchNg21.h>
#include <search_schemes/generator/optimum.h>
#include <search_schemes/expand.h>

namespace search {

search_schemes::Scheme const& search_scheme_cache::get(
    size_t const pex_leaf_query_length,
    size_t const pex_leaf_num_errors
) {
    auto const scheme_data = std::make_tuple(pex_leaf_query_length, pex_leaf_num_errors);
    auto const iter = schemes.find(scheme_data);

    if (iter == schemes.end()) {
        auto search_scheme = search_schemes::expand(
            search_schemes::generator::optimum(0, pex_leaf_num_errors), 
            pex_leaf_query_length
        );

        auto const [emplaced_iter, _] = schemes.emplace(scheme_data, std::move(search_scheme));
        return emplaced_iter->second;
    }

    return iter->second;
}


bool anchor::is_better_than(anchor const& other) {
    size_t const position_difference = position < other.position ?
        other.position - position : position - other.position;

    return num_errors <= other.num_errors && 
        position_difference <= other.num_errors - num_errors;
}

void anchor::mark_for_erasure() {
    num_errors = internal::erase_marker;
}

bool anchor::should_be_erased() const {
    return num_errors == internal::erase_marker;
}

anchors_by_seed_and_reference search_seeds(
    std::vector<seed> const& seeds,
    fmindex const& index,
    search_scheme_cache& scheme_cache,
    size_t const num_reference_sequences
) {
    anchors_by_seed_and_reference anchors(
        seeds.size(), 
        std::vector<std::vector<anchor>>(num_reference_sequences)
    );

    auto const seeds_span = std::span(seeds);

    for (size_t seed_id = 0; seed_id < seeds.size(); ++seed_id) {
        auto const& seed = seeds[seed_id];
        auto const& search_scheme = scheme_cache.get(
            seed.sequence.size(),
            seed.num_errors
        );

        // wrapper for the search interface that expects a range
        auto const seed_single_span = seeds_span.subspan(seed_id, 1) 
            | std::views::transform(&search::seed::sequence);

        auto& anchors_of_seed = anchors[seed_id];

        // if the preorder function inside fmindex search occurs in any profile, we can optimize it away
        fmindex_collection::search_ng21::search(
            index,
            seed_single_span,
            search_scheme,
            [&index, &anchors_of_seed] ([[maybe_unused]] size_t const seed_id, auto cursor, size_t const errors) {                                
                for (auto anchor : cursor) {
                    auto const [reference_id, position] = index.locate(anchor);
                    anchors_of_seed[reference_id].emplace_back(position, errors);
                }
            }
        );
    }

    internal::erase_useless_anchors(anchors);

    return anchors;
}

namespace internal {

void erase_useless_anchors(anchors_by_seed_and_reference & anchors) {
    for (auto & anchors_of_seed : anchors) {
        for (auto & anchors_of_seed_and_reference : anchors_of_seed) {
            // this must stay, otherwise the expression in the below for loop head could underflow
            if (anchors_of_seed_and_reference.empty()) {
                continue;
            }

            std::ranges::sort(anchors_of_seed_and_reference, {}, [] (anchor const& h) { return h.position; });

            for (size_t current_anchor_index = 0; current_anchor_index < anchors_of_seed_and_reference.size() - 1;) {
                auto & current_anchor = anchors_of_seed_and_reference[current_anchor_index];
                size_t other_anchor_index = current_anchor_index + 1;

                while (
                    other_anchor_index < anchors_of_seed_and_reference.size() &&
                    current_anchor.is_better_than(anchors_of_seed_and_reference[other_anchor_index]))
                {
                    anchors_of_seed_and_reference[other_anchor_index].mark_for_erasure();
                    ++other_anchor_index;
                }

                if (
                    other_anchor_index < anchors_of_seed_and_reference.size() &&
                    anchors_of_seed_and_reference[other_anchor_index].is_better_than(current_anchor)) {
                    current_anchor.mark_for_erasure();
                }

                current_anchor_index = other_anchor_index;
            }

            std::erase_if(anchors_of_seed_and_reference, [] (anchor const& a) { return a.should_be_erased(); } );
        }
    }
}

} // namespace internal

} // namespace search
