#pragma once

#include <fmindex.hpp>
#include <tuple_hash.hpp>

#include <string_view>
#include <span>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <search_schemes/Scheme.h>

namespace search {

struct seed {
    std::span<const uint8_t> const sequence;
    size_t const num_errors;
    size_t const query_position;
};

// group count and score are not compared
bool operator==(seed const& lhs, seed const& rhs);

struct anchor_t {
    size_t reference_position;
    size_t num_errors;
    // size_t query_position;
    // size_t length;

    bool is_better_than(anchor_t const& other);

    void mark_for_erasure();

    bool should_be_erased() const;
};

bool operator==(anchor_t const& lhs, anchor_t const& rhs);

using anchors = std::vector<anchor_t>;

enum class anchor_group_order_t {
    num_errors_first, count_first, hybrid
};

anchor_group_order_t anchor_group_order_from_string(std::string_view const s);

struct search_config {
    size_t const max_num_anchors;
    anchor_group_order_t const anchor_group_order;
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

// not thread safe due to search scheme cache
struct searcher {
    fmindex const& index;
    size_t const num_reference_sequences;
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

namespace internal {

struct anchor_group {
    fmindex_cursor cursor;
    size_t num_errors;
};

static inline constexpr size_t erase_marker = std::numeric_limits<size_t>::max();

// returns the number of kept anchors, sorts anchors by position
size_t erase_useless_anchors(std::vector<anchors>& anchors_by_reference);

} // namespace internal

} // namespace search
