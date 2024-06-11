#pragma once

#include <fmindex.hpp>
#include <tuple_hash.hpp>

#include <span>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <search_schemes/Scheme.h>

namespace search {

struct seed {
    std::span<const uint8_t> const sequence;
    size_t const num_errors;
};

struct anchor_t {
    size_t position;
    size_t num_errors;

    bool is_better_than(anchor_t const& other);

    void mark_for_erasure();

    bool should_be_erased() const;
};

using anchors = std::vector<anchor_t>;

struct search_config {
    // if the number of anchors for a seed exceeds this threshold,
    // no anchors will be reported for that seed
    // current problem: the number of raw anchors is larger than what
    // we actually need, because of the repetitive alignments,
    // so currently this is difficult to set right
    size_t const max_num_raw_anchors;

    size_t const max_num_errors;
};

enum class seed_status {
    fully_excluded, partly_excluded, not_excluded
};

struct search_result {
    struct anchors_of_seed {
        seed_status const status;
        size_t const num_kept_useful_anchors;
        size_t const num_excluded_raw_anchors;

        // empty if fully excluded
        std::vector<anchors> const anchors_by_reference;
    };

    std::vector<anchors_of_seed> const anchors_by_seed{};
    size_t const num_fully_excluded_seeds;
};

class search_scheme_cache;

// not thread safe due to search scheme cache
struct searcher {
    fmindex const& index;
    size_t const num_reference_sequences;
    search_scheme_cache& scheme_cache;
    search_config const config;

    search_result search_seeds(
        std::vector<seed> const& seeds
    ) const;
};

class search_scheme_cache {
public:
    search_schemes::Scheme const& get(
        size_t const pex_leaf_query_length,
        size_t const pex_leaf_num_errors
    );

private:
    std::unordered_map<std::tuple<size_t, size_t>, search_schemes::Scheme> schemes;
};

struct cursors_t {
    struct cursors_with_given_num_errors_t {
        std::vector<fmindex_cursor> cursors;
        size_t total_num_raw_anchors{};
    };

    cursors_t(search_config const& config);

    std::vector<cursors_with_given_num_errors_t> cursors_by_num_errors;
    size_t total_num_raw_anchors{};
};

namespace internal {

static inline constexpr size_t erase_marker = std::numeric_limits<size_t>::max();

void erase_useless_anchors(anchors& anchors_of_seed_and_reference);

} // namespace internal

} // namespace search
