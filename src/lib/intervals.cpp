#include <intervals.hpp>
#include <math.hpp>

#include <algorithm>
#include <cassert>
#include <mutex>
#include <stdexcept>
#include <vector>

namespace intervals {

size_t half_open_interval::size() const {
    return end - start;
}

half_open_interval half_open_interval::overlap_interval_with(half_open_interval const other) const {
    assert(relationship_with(other) != interval_relationship::completely_above);
    assert(relationship_with(other) != interval_relationship::completely_below);

    return half_open_interval {
        .start = std::max(start, other.start),
        .end = std::min(end, other.end)
    };
}

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

lib_interval_tree::interval<size_t> half_open_interval::to_lib_intervaltree_interval() const {
    return { start, end };
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

void verified_intervals::configure(
    use_interval_optimization const activity_status_,
    double const overlap_rate_that_counts_as_contained_
) {
    activity_status = activity_status_;
    overlap_rate_that_counts_as_contained = overlap_rate_that_counts_as_contained_;
}

void verified_intervals::insert(half_open_interval const new_interval) {
    if (activity_status == use_interval_optimization::off) {
        return;
    }

    intervals.insert_overlap(new_interval.to_lib_intervaltree_interval());
    intervals.deoverlap();
}

bool verified_intervals::contains(
    half_open_interval const target_interval
) const {
    if (activity_status == use_interval_optimization::off) {
        return false;
    }

    bool does_contain = false;

    intervals.overlap_find_all(target_interval.to_lib_intervaltree_interval(), [&] (auto interval_it) {
        half_open_interval const existing_interval = {
            .start = interval_it->low(),
            .end = interval_it->high()
        };

        // the return value indicates whether to keep searching
        switch (existing_interval.relationship_with(target_interval)) {
            case interval_relationship::completely_above:
            case interval_relationship::completely_below:
                throw std::runtime_error("(should not happen) found disjoint interval in overlap function");
            case interval_relationship::equal:
            case interval_relationship::contains:
                does_contain = true;
                break;
            case interval_relationship::inside:
            case interval_relationship::overlapping_or_touching_above:
            case interval_relationship::overlapping_or_touching_below:
                static constexpr double epsilon = 0.000000001;
                does_contain = static_cast<double>(
                    target_interval.overlap_interval_with(existing_interval).size()
                ) / target_interval.size() + epsilon >= overlap_rate_that_counts_as_contained;
        }

        return !does_contain;
    });

    return does_contain;
}

std::shared_ptr<verified_intervals_for_all_references> create_thread_safe_verified_intervals(
    size_t const num_references,
    use_interval_optimization const activity_status,
    double const overlap_rate_that_counts_as_contained
) {
    auto out = std::make_shared<verified_intervals_for_all_references>(num_references);

    for (size_t i = 0; i < num_references; ++i) {
        auto && [lock, ivls] = out->at(i).lock_unique();
        ivls.configure(activity_status, overlap_rate_that_counts_as_contained);
    }

    return out;
}

} // namespace intervals
