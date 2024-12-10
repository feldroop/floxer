#pragma once

#include <alignment.hpp>
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
    size_t const pex_leaf_index;
};

// group count and score are not compared
bool operator==(seed const& lhs, seed const& rhs);

struct anchor_t {
    size_t pex_leaf_index; // a.k.a. seed id
    size_t reference_id;
    size_t reference_position;
    size_t num_errors;

    bool is_better_than(anchor_t const& other);

    void mark_for_erasure();

    bool should_be_erased() const;
};

bool operator==(anchor_t const& lhs, anchor_t const& rhs);

using anchors_t = std::vector<anchor_t>;

enum class anchor_group_order_t {
    num_errors_first, count_first, none
};

anchor_group_order_t anchor_group_order_from_string(std::string_view const s);

enum class anchor_choice_strategy_t {
    round_robin, full_groups
};

anchor_choice_strategy_t anchor_choice_strategy_from_string(std::string_view const s);

struct search_config {
    size_t const max_num_anchors_hard;
    size_t const max_num_anchors_soft;
    anchor_group_order_t const anchor_group_order;
    anchor_choice_strategy_t const anchor_choice_strategy;
};

struct anchor_package {
    size_t package_id;
    anchors_t anchors;
    alignment::query_orientation orientation;
};

struct search_result {
    struct anchor_iterator {
        search_result const& res;
        size_t current_seed_index;
        size_t current_reference_index;
        size_t current_anchor_index;

        std::optional<std::reference_wrapper<const anchor_t>> next();
    };

    struct anchors_of_seed {
        size_t num_kept_useful_anchors;
        size_t num_kept_raw_anchors;
        size_t num_excluded_raw_anchors_by_soft_cap;

        // empty if fully excluded
        std::vector<anchors_t> anchors_by_reference;
    };

    std::vector<anchors_of_seed> anchors_by_seed{};
    size_t num_fully_excluded_seeds;

    // a "flattened" iterator over the anchors of all queries to all references
    // the order is sorted by query first, then by reference, then by position (the only way it makes sense)
    anchor_iterator anchor_iter() const;

    // package for verification tasks
    void append_anchor_packages(
        std::vector<anchor_package>& out_packages,
        size_t const num_anchors_per_package,
        alignment::query_orientation orientation
    ) const;
};

struct searcher {
    fmindex const& index;
    size_t const num_reference_sequences;
    search_config const config;

    search_result search_seeds(
        std::vector<seed> const& seeds
    ) const;
};

namespace internal {

class search_scheme_cache {
public:
    search_schemes::Scheme const& get(
        size_t const pex_leaf_query_length,
        size_t const pex_leaf_num_errors
    );

private:
    std::unordered_map<std::tuple<size_t, size_t>, search_schemes::Scheme> schemes;
};

struct anchor_group {
    fmindex_cursor cursor;
    size_t num_errors;
};

static inline constexpr size_t erase_marker = std::numeric_limits<size_t>::max();

// returns the number of kept anchors, sorts anchors by position
size_t erase_useless_anchors(std::vector<anchors_t>& anchors_by_reference);

} // namespace internal

} // namespace search
