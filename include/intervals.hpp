#pragma once

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

// important invariant: the intervals stored are disjoint
class interval_set {
public:
    void insert(half_open_interval const new_interval);

    bool contains(half_open_interval const target_interval) const;

    size_t size() const;

private:
    std::set<half_open_interval> intervals{};
};

} // namespace intervals
