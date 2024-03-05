#include <fmindex.hpp>
#include <pex.hpp>
#include <search.hpp>

#include <algorithm>
#include <limits>

#include <fmindex-collection/search/SearchNg21.h>
#include <search_schemes/generator/optimum.h>
#include <search_schemes/expand.h>

namespace search {

search_scheme_cache::search_scheme_cache(size_t const pex_leaf_num_errors_)
        : pex_leaf_num_errors{pex_leaf_num_errors_} {}

search_schemes::Scheme const& search_scheme_cache::get(size_t const pex_leaf_query_length) {
    if (!schemes.contains(pex_leaf_query_length)) {
        auto search_scheme = search_schemes::expand(
            search_schemes::generator::optimum(0, pex_leaf_num_errors), 
            pex_leaf_query_length
        );

        schemes.emplace(pex_leaf_query_length, std::move(search_scheme));
    }

    return schemes.at(pex_leaf_query_length);
}


bool hit::is_better_than(hit const& other) {
    size_t const position_difference = position < other.position ?
        other.position - position : position - other.position;

    return num_errors <= other.num_errors && 
        position_difference <= other.num_errors - num_errors;
}

void erase_useless_hits(hit_list & hits) {
    static constexpr size_t erase_marker = std::numeric_limits<size_t>::max();
    
    for (auto & query_hits : hits) {
        for (auto & query_to_ref_hits : query_hits) {
            // this must stay, otherwise the expression in the below for loop head could underflow
            if (query_to_ref_hits.empty()) {
                continue;
            }

            std::ranges::sort(query_to_ref_hits, {}, [] (hit const& h) { return h.position; });

            for (size_t i = 0; i < query_to_ref_hits.size() - 1; ++i) {
                auto & current_hit = query_to_ref_hits[i];
                size_t j = i + 1;

                while (j < query_to_ref_hits.size() && current_hit.is_better_than(query_to_ref_hits[j])) { 
                    query_to_ref_hits[j].num_errors = erase_marker;
                    ++j; 
                }

                if (j < query_to_ref_hits.size() && query_to_ref_hits[j].is_better_than(current_hit)) {
                    current_hit.num_errors = erase_marker;
                }
            }

            std::erase_if(query_to_ref_hits, [] (hit const& h) { return h.num_errors == erase_marker; } );
        }
    }
}

hit_list search_leaf_queries(
    std::vector<std::span<const uint8_t>> const& leaf_queries,
    fmindex const& index,
    search_schemes::Scheme const& search_scheme,
    size_t const num_reference_sequences
) {    
    hit_list hits(
        leaf_queries.size(), 
        std::vector<std::vector<hit>>(num_reference_sequences)
    );

    fmindex_collection::search_ng21::search(
        index,
        leaf_queries,
        search_scheme,
        [&index, &hits] (size_t const leaf_query_id, auto cursor, size_t const errors) {                
            auto& query_hits = hits[leaf_query_id];

            for (auto hit{begin(cursor)}; hit < end(cursor); ++hit) {
                auto const [reference_id, pos] = index.locate(hit);
                query_hits[reference_id].emplace_back(pos, errors);
            }
        }
    );

    erase_useless_hits(hits);

    return hits;
}

} // namespace search
