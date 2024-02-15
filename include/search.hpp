#pragma once 

#include <vector>
#include <unordered_map>

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
};

// REFACTOR LATER hits[leaf_query_id][reference_id] -> hits
using hit_list = std::vector<std::vector<std::vector<hit>>>;

hit_list search_fastq_query(
    std::vector<uint8_t> const& fastq_query,
    fmindex const& index,
    pex_tree const& tree,
    search_scheme_cache& scheme_cache,
    size_t const num_reference_sequences
);

} // namespace search
