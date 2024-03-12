#pragma once 

#include <miscellaneous.hpp>

#include <span>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <search_schemes/Scheme.h>

namespace search {

class search_scheme_cache {
public:    
    search_schemes::Scheme const& get(
        size_t const pex_leaf_query_length,
        size_t const pex_leaf_num_errors
    );

private:
    std::unordered_map<std::tuple<size_t, size_t>, search_schemes::Scheme> schemes;
};

struct query {
    std::span<const uint8_t> const sequence;
    size_t const num_errors;
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
    std::vector<query> const& leaf_queries,
    fmindex const& index,
    search_scheme_cache& scheme_cache,
    size_t const num_reference_sequences
);

} // namespace search
