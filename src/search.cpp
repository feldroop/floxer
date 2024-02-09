#include <search.hpp>

#include <search_schemes/generator/optimum.h>
#include <search_schemes/expand.h>

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
