#pragma once

#include <alignment.hpp>

#include <cstddef>
#include <optional>
#include <set>

namespace intervals {

enum interval_relationship {
    completely_above,
    completely_below,
    contains,
    equal,
    inside,
    overlapping_or_touching_above,
    overlapping_or_touching_below
};

// [start, end), must be non-empty (start < end)
struct half_open_interval {
    size_t start;
    size_t end;

    // <this> has returned relationship with <other>
    interval_relationship relationship_with(half_open_interval const other) const;
};

// order intervals by their END position, because it works better with the std::set functions
auto operator<=>(half_open_interval const& ivl1, half_open_interval const& ivl2);

// important invariants: both of the interval sets store only disjoint intervals
class verified_intervals {
public:
    using intervals_t = std::set<half_open_interval>;

    void insert(half_open_interval const new_interval, alignment::alignment_outcome const outcome);

    // std::nullopt if it's not contained
    std::optional<alignment::alignment_outcome> contains(half_open_interval const target_interval) const;

    size_t size() const;

private:
    void given_set_insert(half_open_interval const new_interval, intervals_t& intervals);

    bool given_set_contains(half_open_interval const target_interval, intervals_t const& intervals) const;

    intervals_t intervals_with_alignment{};
    intervals_t intervals_without_alignment{};
};

} // namespace intervals
