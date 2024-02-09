#pragma once 

#include <search_schemes/Scheme.h>

#include <unordered_map>

class search_scheme_cache {
public:
    search_scheme_cache() = delete;
    search_scheme_cache(size_t const pex_leaf_num_errors_);
    
    search_schemes::Scheme const& get(size_t const query_length);
private:
    size_t const pex_leaf_num_errors;

    // query length determines search scheme structure uniquely for now
    // actually it also depends on the number of errors
    std::unordered_map<size_t, search_schemes::Scheme> schemes;
};
