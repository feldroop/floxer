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

    half_open_interval trim_from_both_sides(size_t const amount) const;
};

bool operator==(half_open_interval const& interval1, half_open_interval const& interval2);

// order intervals by their END position, because it works better with the std::set functions
auto operator<(half_open_interval const& interval1, half_open_interval const& interval2);

enum use_interval_optimization {
    on, off
};

// important invariants: this set stores only disjoint intervals
class verified_intervals {
public:
    verified_intervals(use_interval_optimization const activity_status);

    using intervals_t = std::set<half_open_interval>;

    void insert(half_open_interval const new_interval);

    bool contains(half_open_interval const target_interval) const;

    size_t size() const;

private:
    use_interval_optimization const activity_status;

    std::set<half_open_interval> intervals{};
};

} // namespace intervals
