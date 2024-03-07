#pragma once 

#include <span>
#include <unordered_map>
#include <vector>

#include <search_schemes/Scheme.h>

namespace search {

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

struct hit {
    size_t position;
    size_t num_errors;

    bool is_better_than(hit const& other);
    void mark_for_erasure();
    bool should_be_erased() const;
};

// REFACTOR LATER hits[leaf_query_id][reference_id] -> hits
using hit_list = std::vector<std::vector<std::vector<hit>>>;

hit_list search_leaf_queries(
    std::vector<std::span<const uint8_t>> const& leaf_queries,
    fmindex const& index,
    search_schemes::Scheme const& search_scheme,
    size_t const num_reference_sequences
);

} // namespace search
