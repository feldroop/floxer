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


bool hit::is_better_than(hit const& other) {
    size_t const position_difference = position < other.position ?
        other.position - position : position - other.position;

    return num_errors <= other.num_errors && 
        position_difference <= other.num_errors - num_errors;
}

static constexpr size_t erase_marker = std::numeric_limits<size_t>::max();

void hit::mark_for_erasure() {
    num_errors = erase_marker;
}

bool hit::should_be_erased() const {
    return num_errors == erase_marker;
}

void erase_useless_hits(hit_list & hits) {
    for (auto & query_hits : hits) {
        for (auto & query_to_ref_hits : query_hits) {
            // this must stay, otherwise the expression in the below for loop head could underflow
            if (query_to_ref_hits.empty()) {
                continue;
            }

            std::ranges::sort(query_to_ref_hits, {}, [] (hit const& h) { return h.position; });

            for (size_t i = 0; i < query_to_ref_hits.size() - 1;) {
                auto & current_hit = query_to_ref_hits[i];
                size_t j = i + 1;

                while (j < query_to_ref_hits.size() && current_hit.is_better_than(query_to_ref_hits[j])) {
                    query_to_ref_hits[j].mark_for_erasure();
                    ++j;
                }

                if (j < query_to_ref_hits.size() && query_to_ref_hits[j].is_better_than(current_hit)) {
                    current_hit.mark_for_erasure();
                }

                i = j;
            }

            std::erase_if(query_to_ref_hits, [] (hit const& h) { return h.should_be_erased(); } );
        }
    }
}

hit_list search_leaf_queries(
    std::vector<query> const& leaf_queries,
    fmindex const& index,
    search_scheme_cache& scheme_cache,
    size_t const num_reference_sequences
) {
    hit_list hits(
        leaf_queries.size(), 
        std::vector<std::vector<hit>>(num_reference_sequences)
    );

    auto const leaf_queries_span = std::span(leaf_queries);

    for (size_t i = 0; i < leaf_queries.size(); ++i) {
        auto const& leaf_query = leaf_queries[i];
        auto const& search_scheme = scheme_cache.get(
            leaf_query.sequence.size(),
            leaf_query.num_errors
        );

        // wrapper for the search interface that expects something range-like
        auto const leaf_query_single_span = leaf_queries_span.subspan(i, 1) 
            | std::views::transform(&search::query::sequence);

        // if the preorder function inside fmindex search occurs in any perf profile, we can optimize it
        // also, search_ng22 is for the aidditional alignment output
        fmindex_collection::search_ng21::search(
            index,
            leaf_query_single_span,
            search_scheme,
            [&index, &hits] (size_t const leaf_query_id, auto cursor, size_t const errors) {                
                auto& query_hits = hits[leaf_query_id];

                for (auto hit{begin(cursor)}; hit < end(cursor); ++hit) {
                    auto const [reference_id, pos] = index.locate(hit);
                    query_hits[reference_id].emplace_back(pos, errors);
                }
            }
        );
    }

    erase_useless_hits(hits);

    return hits;
}

} // namespace search
