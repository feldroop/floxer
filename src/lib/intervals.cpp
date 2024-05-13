#include <intervals.hpp>

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <vector>

namespace intervals {

interval_relationship half_open_interval::relationship_with(half_open_interval const other) const {
    assert(start < end);
    assert(other.start < other.end);

    if (start > other.end) {
        return interval_relationship::completely_above;
    } else if (end < other.start) {
        return interval_relationship::completely_below;
    } else if (start == other.start && end == other.end) {
        return interval_relationship::equal;
    } else if (start <= other.start && end >= other.end) {
        return interval_relationship::contains;
    } else if (start >= other.start && end <= other.end) {
        return interval_relationship::inside;
    } else if (start > other.start && start <= other.end) {
        return interval_relationship::overlapping_or_touching_above;
    } else {
        assert(end < other.end && end >= other.start);
        return interval_relationship::overlapping_or_touching_below;
    }
}

half_open_interval half_open_interval::trim_from_both_sides(size_t const amount) const {
    assert(start < end);

    size_t const new_end = std::max(start + 1, amount > end ? 0 : end - amount);
    size_t const new_start = std::min(new_end - 1, start + amount);

    return half_open_interval {
        .start = new_start,
        .end = new_end
    };
}

bool operator==(half_open_interval const& interval1, half_open_interval const& interval2) {
    assert(interval1.start < interval1.end);
    assert(interval2.start < interval2.end);

    return interval1.start == interval2.start && interval1.end == interval2.end;
}

auto operator<(half_open_interval const& interval1, half_open_interval const& interval2) {
    assert(interval1.start < interval1.end);
    assert(interval2.start < interval2.end);

    return interval1.end < interval2.end;
}

verified_intervals::verified_intervals(use_interval_optimization const activity_status_)
    : activity_status{activity_status_} {}

void verified_intervals::insert(half_open_interval const new_interval, alignment::alignment_outcome const outcome) {
    if (activity_status == use_interval_optimization::off) {
        return;
    }

    switch (outcome) {
        case alignment::alignment_outcome::alignment_exists:
            given_set_insert(new_interval, intervals_with_alignment);
            return;
        case alignment::alignment_outcome::no_adequate_alignment_exists:
            given_set_insert(new_interval, intervals_without_alignment);
            return;
        default:
            throw std::runtime_error("(should be unreachable) internal bug in interval set");
    }
}

std::optional<alignment::alignment_outcome> verified_intervals::contains(half_open_interval const target_interval) const {
    if (activity_status == use_interval_optimization::off) {
        return std::nullopt;
    }

    // it is important that we ask first for the alignment_exists case, because it should override the other case
    if (given_set_contains(target_interval, intervals_with_alignment)) {
        return alignment::alignment_outcome::alignment_exists;
    } else if (given_set_contains(target_interval, intervals_without_alignment)) {
        return alignment::alignment_outcome::no_adequate_alignment_exists;
    } else {
        return std::nullopt;
    }
}

void verified_intervals::given_set_insert(half_open_interval const new_interval, verified_intervals::intervals_t& intervals) {
    if (intervals.empty()) {
        intervals.insert(new_interval);
        return;
    }

    std::vector<half_open_interval> intervals_to_remove{};
    half_open_interval interval_to_insert = new_interval;

    // first interval that has an end position ABOVE new_interval.end
    auto existing_interval_iter = intervals.upper_bound(new_interval);

    if (existing_interval_iter == intervals.end()) {
        // safe because intervals can't be empty
        --existing_interval_iter;
    }

    bool continue_searching = true;

    // This is linear time in the worst case, if we have to merge all intervals.
    // However, I believe this case will rarely happen in the usage of this program.
    // Also it allows simple and quick contains queries, which is a nice trade-off.
    while (continue_searching) {
        auto const existing_interval = *existing_interval_iter;
        auto const relationship = existing_interval.relationship_with(new_interval);

        switch (relationship) {
            case interval_relationship::completely_above:
                break;
            case interval_relationship::completely_below:
                continue_searching = false;
                break;
            case interval_relationship::contains:
            case interval_relationship::equal:
                return;
            case interval_relationship::inside:
                intervals_to_remove.emplace_back(existing_interval);
                break;
            case interval_relationship::overlapping_or_touching_above:
                intervals_to_remove.emplace_back(existing_interval);
                interval_to_insert.end = existing_interval.end;
                break;
            case interval_relationship::overlapping_or_touching_below:
                intervals_to_remove.emplace_back(existing_interval);
                interval_to_insert.start = existing_interval.start;
                continue_searching = false;
                break;
            default:
                throw std::runtime_error("(should be unreachable) internal bug in interval set");
        }

        if (existing_interval_iter == intervals.begin()) {
            continue_searching = false;
        } else {
            --existing_interval_iter;
        }
    }

    for (half_open_interval const interval_to_remove : intervals_to_remove) {
        [[maybe_unused]] size_t const num_erased = intervals.erase(interval_to_remove);
        assert(num_erased == 1);
    }

    intervals.insert(interval_to_insert);
}

bool verified_intervals::given_set_contains(
    half_open_interval const target_interval,
    verified_intervals::intervals_t const& intervals
) const {
    // first interval that has an end position NOT BELOW target_interval.end
    auto const existing_interval_iter = intervals.lower_bound(target_interval);

    if (existing_interval_iter == intervals.end()) {
        // all existing intervals have an end position below target interval
        return false;
    }

    auto const relationship = existing_interval_iter->relationship_with(target_interval);

    return relationship == interval_relationship::contains || relationship == interval_relationship::equal;
}

size_t verified_intervals::size() const {
    return intervals_with_alignment.size() + intervals_without_alignment.size();
}

} // namespace intervals
