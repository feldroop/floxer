#include <fmindex.hpp>
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

search_result searcher::search_seeds(
    std::vector<seed> const& seeds
) const {
    std::vector<search_result::anchors_of_seed> anchors_by_seed{};
    size_t num_excluded_seeds = 0;

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

        size_t num_raw_anchors = 0;
        std::vector<internal::fmindex_search_return> fmindex_search_returns{};

        // if the preorder function inside fmindex search occurs in any profile, we can optimize it away
        fmindex_collection::search_ng21::search(
            index,
            seed_single_span,
            search_scheme,
            [&fmindex_search_returns, &num_raw_anchors]
            ([[maybe_unused]] size_t const seed_id, auto cursor, size_t const errors) {
                num_raw_anchors += cursor.count();
                fmindex_search_returns.emplace_back(internal::fmindex_search_return{
                    .cursor = std::move(cursor),
                    .num_errors = errors
                });
                // as a possible hack: could throw here to stop the search as future optimization
            }
        );

        if (num_raw_anchors >= config.max_num_raw_anchors) {
            anchors_by_seed.emplace_back(search_result::anchors_of_seed {
                .excluded = true,
                .num_anchors = num_raw_anchors,
                .anchors_by_reference = std::vector<anchors>{}
            });
            continue;
        }

        std::vector<anchors> anchors_by_reference(num_reference_sequences);
        for (auto const& [cursor, num_errors] : fmindex_search_returns) {
            for (auto const& anchor: cursor) {
                auto const [reference_id, position] = index.locate(anchor);
                anchors_by_reference[reference_id].emplace_back(position, num_errors);
            }
        }

        size_t num_useful_anchors = 0;
        for (auto& anchors_of_seed_and_reference : anchors_by_reference) {
            internal::erase_useless_anchors(anchors_of_seed_and_reference);
            num_useful_anchors += anchors_of_seed_and_reference.size();
        }

        anchors_by_seed.emplace_back(search_result::anchors_of_seed{
            .excluded = false,
            .num_anchors = num_useful_anchors,
            .anchors_by_reference = std::move(anchors_by_reference)
        });
    }

    return search_result {
        .anchors_by_seed = std::move(anchors_by_seed),
        .num_excluded_seeds = num_excluded_seeds
    };
}

namespace internal {

void erase_useless_anchors(anchors& anchors_of_seed_and_reference) {
    // this must stay, otherwise the expression in the below for loop head could underflow
    if (anchors_of_seed_and_reference.empty()) {
        return;
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

} // namespace internal

} // namespace search
