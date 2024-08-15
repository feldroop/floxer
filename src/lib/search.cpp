#include <fmindex.hpp>
#include <search.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <ranges>

#include <fmindex-collection/search/SearchNg21.h>
#include <search_schemes/generator/optimum.h>
#include <search_schemes/generator/h2.h>
#include <search_schemes/expand.h>

namespace search {

bool operator==(seed const& lhs, seed const& rhs) {
    if (lhs.sequence.size() != rhs.sequence.size()) {
        return false;
    }

    for (size_t i = 0; i < lhs.sequence.size(); ++i) {
        if (lhs.sequence[i] != rhs.sequence[i]) {
            return false;
        }
    }

    return lhs.num_errors == rhs.num_errors &&
        lhs.query_position == rhs.query_position;
}

search_schemes::Scheme const& search_scheme_cache::get(
    size_t const pex_leaf_query_length,
    size_t const pex_leaf_num_errors
) {
    auto const scheme_data = std::make_tuple(pex_leaf_query_length, pex_leaf_num_errors);
    auto const iter = schemes.find(scheme_data);

    if (iter == schemes.end()) {
        auto search_scheme = search_schemes::expand(
            (pex_leaf_num_errors <= 3) ?
                search_schemes::generator::optimum(0, pex_leaf_num_errors) :
                // h2 = heuristic 2, the best heuristic search scheme generator,
                // because the optima are not known for more than 3 errors
                search_schemes::generator::h2(pex_leaf_num_errors + 2, 0, pex_leaf_num_errors),
            pex_leaf_query_length
        );

        auto const [emplaced_iter, _] = schemes.emplace(scheme_data, std::move(search_scheme));
        return emplaced_iter->second;
    }

    return iter->second;
}

bool anchor_t::is_better_than(anchor_t const& other) {
    size_t const position_difference = reference_position < other.reference_position ?
        other.reference_position - reference_position : reference_position - other.reference_position;

    return num_errors <= other.num_errors &&
        position_difference <= other.num_errors - num_errors;
}

void anchor_t::mark_for_erasure() {
    num_errors = internal::erase_marker;
}

bool anchor_t::should_be_erased() const {
    return num_errors == internal::erase_marker;
}

anchor_group_order_t anchor_group_order_from_string(std::string_view const s) {
    if (s == "errors_first") {
        return anchor_group_order_t::num_errors_first;
    } else if (s == "count_first") {
        return anchor_group_order_t::count_first;
    } else {
        return anchor_group_order_t::hybrid;
    }
}

search_result searcher::search_seeds(
    std::vector<seed> const& seeds
) const {
    using namespace internal;

    std::vector<search_result::anchors_of_seed> anchors_by_seed{};
    size_t num_fully_excluded_seeds = 0;

    auto const seeds_span = std::span(seeds);

    for (size_t seed_id = 0; seed_id < seeds.size(); ++seed_id) {
        auto const& seed = seeds[seed_id];
        auto const& search_scheme = scheme_cache.get(
            seed.sequence.size(),
            seed.num_errors
        );

        // wrapper for the search interface that expects a range
        auto const seed_single_span = seeds_span.subspan(seed_id, 1)
            | std::views::transform(&search::seed::sequence);

        std::vector<anchor_group> anchor_groups{};
        size_t total_num_raw_anchors = 0;

        fmindex_collection::search_ng21::search(
            index,
            seed_single_span,
            search_scheme,
            [&anchor_groups, &total_num_raw_anchors] ([[maybe_unused]] size_t const seed_id, auto cursor, size_t const errors) {
                anchor_groups.emplace_back(cursor, errors);
                total_num_raw_anchors += cursor.count();
            }
        );

        switch (config.anchor_group_order) {
            case anchor_group_order_t::count_first:
                std::ranges::sort(anchor_groups, [] (anchor_group const& group1, anchor_group const& group2) {
                    if (group1.cursor.count() != group2.cursor.count()) {
                        return group1.cursor.count() < group2.cursor.count();
                    } else {
                        return group1.num_errors < group2.num_errors;
                    }
                });
                break;

            case anchor_group_order_t::num_errors_first:
                std::ranges::sort(anchor_groups, [] (anchor_group const& group1, anchor_group const& group2) {
                    if (group1.num_errors != group2.num_errors) {
                        return group1.cursor.count() < group2.cursor.count();
                    } else {
                        return group1.num_errors < group2.num_errors;
                    }
                });
                break;

            case anchor_group_order_t::hybrid:
                std::ranges::sort(anchor_groups, [] (anchor_group const& group1, anchor_group const& group2) {
                    return (group1.num_errors + 1) * group1.cursor.count() <
                        (group2.num_errors + 1) * group2.cursor.count();
                });
                break;

            default:
                throw std::runtime_error("(Should be unreachable) internal bug in anchor group order config.");
        }

        size_t num_kept_raw_anchors = 0;
        std::vector<anchors> anchors_by_reference(num_reference_sequences);
        size_t i = 0;
        while (
            i < anchor_groups.size() &&
            num_kept_raw_anchors + anchor_groups[i].cursor.count() <= config.max_num_located_raw_anchors
        ) {
            auto const& cursor = anchor_groups[i].cursor;
            size_t const group_count = cursor.count();

            for (auto const& anchor: cursor) {
                auto const [reference_id, position] = index.locate(anchor);
                anchors_by_reference[reference_id].emplace_back(anchor_t {
                    .reference_position = position,
                    .num_errors = anchor_groups[i].num_errors,
                    .query_position = seed.query_position,
                    .length = seed.sequence.size(),
                    .group_count = group_count
                });
            }

            num_kept_raw_anchors += group_count;
            ++i;
        }

        size_t const num_excluded_raw_anchors = total_num_raw_anchors - num_kept_raw_anchors;

        size_t num_kept_useful_anchors = erase_useless_anchors(anchors_by_reference);

        if (num_kept_useful_anchors > config.max_num_kept_anchors) {
            num_kept_useful_anchors = erase_low_scoring_anchors(anchors_by_reference, config.max_num_kept_anchors);
        }

        seed_status status;
        if (num_kept_useful_anchors == 0) {
            status = seed_status::fully_excluded;
        } else if (num_excluded_raw_anchors == 0 && num_kept_useful_anchors == num_kept_raw_anchors) {
            status = seed_status::not_excluded;
        } else {
            status = seed_status::partly_excluded;
        }

        anchors_by_seed.emplace_back(search_result::anchors_of_seed{
            .status = status,
            .num_kept_useful_anchors = num_kept_useful_anchors,
            .num_excluded_raw_anchors = num_excluded_raw_anchors,
            .anchors_by_reference = std::move(anchors_by_reference)
        });
    }

    return search_result {
        .anchors_by_seed = std::move(anchors_by_seed),
        .num_fully_excluded_seeds = num_fully_excluded_seeds
    };
}

namespace internal {

size_t erase_useless_anchors(std::vector<anchors>& anchors_by_reference) {
    size_t num_kept_useful_anchors = 0;

    for (auto& anchors_of_seed_and_reference : anchors_by_reference) {
        // this must stay, otherwise the expression in the below for loop head could underflow
        if (anchors_of_seed_and_reference.empty()) {
            continue;
        }

        std::ranges::sort(anchors_of_seed_and_reference, {}, [] (anchor_t const& h) { return h.reference_position; });

        for (size_t current_anchor_index = 0; current_anchor_index < anchors_of_seed_and_reference.size() - 1;) {
            auto & current_anchor = anchors_of_seed_and_reference[current_anchor_index];
            size_t other_anchor_index = current_anchor_index + 1;

            while (
                other_anchor_index < anchors_of_seed_and_reference.size() &&
                current_anchor.is_better_than(anchors_of_seed_and_reference[other_anchor_index])
            ) {
                anchors_of_seed_and_reference[other_anchor_index].mark_for_erasure();
                ++other_anchor_index;
            }

            if (
                other_anchor_index < anchors_of_seed_and_reference.size() &&
                anchors_of_seed_and_reference[other_anchor_index].is_better_than(current_anchor)) {
                current_anchor.mark_for_erasure();
            }

            current_anchor_index = other_anchor_index;
        }

        std::erase_if(anchors_of_seed_and_reference, [] (anchor_t const& a) { return a.should_be_erased(); } );
        num_kept_useful_anchors += anchors_of_seed_and_reference.size();
    }

    return num_kept_useful_anchors;
}

size_t erase_low_scoring_anchors(std::vector<anchors>& anchors_by_reference, size_t const max_num_kept_anchors) {
    std::vector<double> scores{};

    for (auto& anchors_of_seed_and_reference : anchors_by_reference) {
        size_t const num_anchors_of_reference = anchors_of_seed_and_reference.size();
        for (size_t i = 0; i < num_anchors_of_reference; ++i) {
            auto& anchor = anchors_of_seed_and_reference[i];

            size_t const difference_to_left = i == 0 ?
                std::numeric_limits<size_t>::max() :
                anchor.reference_position - anchors_of_seed_and_reference[i - 1].reference_position;

            size_t const difference_to_right = i == num_anchors_of_reference - 1 ?
                std::numeric_limits<size_t>::max() :
                anchors_of_seed_and_reference[i + 1].reference_position - anchor.reference_position;

            double const difference_to_nearest_neighbor = std::min(difference_to_left, difference_to_right);

            // small score is bad
            double const score = num_anchors_of_reference == 1 ? 0.0 :
                difference_to_nearest_neighbor * (1 / std::sqrt(anchor.group_count)) *
                std::sqrt(num_anchors_of_reference);

            anchor.score = score;
            scores.push_back(score);
        }
    }

    std::ranges::nth_element(scores, scores.begin() + max_num_kept_anchors, std::greater<double>());
    double const score_lower_bound = scores.at(max_num_kept_anchors);

    size_t num_kept_anchors = 0;

    for (auto& anchors_of_seed_and_reference : anchors_by_reference) {
        std::erase_if(anchors_of_seed_and_reference, [score_lower_bound] (anchor_t const& a) {
            return a.score.value() <= score_lower_bound;
        } );
        num_kept_anchors += anchors_of_seed_and_reference.size();
    }

    return num_kept_anchors;
}

} // namespace internal

} // namespace search
