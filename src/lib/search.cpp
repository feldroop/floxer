#include <fmindex.hpp>
#include <search.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <ranges>
#include <set>

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

bool operator==(anchor_t const& lhs, anchor_t const& rhs) {
    return lhs.reference_position == rhs.reference_position &&
        lhs.num_errors == rhs.num_errors;
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
    } else if (s == "none") {
        return anchor_group_order_t::none;
    } else {
        throw std::runtime_error("unexpected anchor group order value");
    }
}

anchor_choice_strategy_t anchor_choice_strategy_from_string(std::string_view const s) {
    if (s == "round_robin") {
        return anchor_choice_strategy_t::round_robin;
    } else if (s == "full_groups") {
        return anchor_choice_strategy_t::full_groups;
    } else if (s == "first_reported") {
        return anchor_choice_strategy_t::first_reported;
    } else {
        throw std::runtime_error("unexpected anchor choice strategy value");
    }
}

std::optional<std::reference_wrapper<const anchor_t>> search_result::anchor_iterator::next() {
    while (current_seed_index < res.anchors_by_seed.size()) {
        auto const& current_anchors_by_seed = res.anchors_by_seed[current_seed_index];

        while (current_reference_index < current_anchors_by_seed.anchors_by_reference.size()) {
            auto const& current_anchors_by_seed_and_reference = current_anchors_by_seed.anchors_by_reference[current_reference_index];

            if (current_anchor_index < current_anchors_by_seed_and_reference.size()) {
                ++current_anchor_index;
                return std::make_optional(std::cref(current_anchors_by_seed_and_reference[current_anchor_index - 1]));
            }

            ++current_reference_index;
            current_anchor_index = 0;
        }

        ++current_seed_index;
        current_reference_index = 0;
        current_anchor_index = 0;
    }

    return std::nullopt;
}

search_result::anchor_iterator search_result::anchor_iter() const {
    return anchor_iterator {
        .res = *this,
        .current_seed_index = 0,
        .current_reference_index = 0,
        .current_anchor_index = 0
    };
}

void search_result::append_anchor_packages(
    std::vector<anchor_package>& out_packages,
    size_t const num_anchors_per_package,
    alignment::query_orientation orientation
) const {
    auto iter = anchor_iter();
    bool anchors_remaining = true;

    while (anchors_remaining) {
        auto package = anchor_package {
            .package_id = out_packages.size(),
            .anchors = anchors_t(),
            .orientation = orientation
        };

        while (package.anchors.size() < num_anchors_per_package) {
            auto anchor_opt = iter.next();

            if(!anchor_opt) {
                anchors_remaining = false;
                break;
            }

            package.anchors.emplace_back(anchor_opt->get());
        }

        if (!package.anchors.empty()) {
            out_packages.emplace_back(std::move(package));
        }
    }
}

search_result searcher::search_seeds(
    std::vector<seed> const& seeds
) const {
    using namespace internal;

    std::vector<search_result::anchors_of_seed> anchors_by_seed{};
    size_t num_fully_excluded_seeds = 0;

    auto const seeds_span = std::span(seeds);
    // this scheme cache exists, because the seeds are not necessarily the same length.
    // the creation of the expanded search scheme in not completely free, therefore we
    // want to reuse them as much as possible. But to be honest, this likely doesn't really
    // matter for the running time of the whole program.
    internal::search_scheme_cache scheme_cache;

    for (size_t seed_index = 0; seed_index < seeds.size(); ++seed_index) {
        auto const& seed = seeds[seed_index];
        auto const& search_scheme = scheme_cache.get(
            seed.sequence.size(),
            seed.num_errors
        );

        // wrapper for the search interface that expects a range
        auto const seed_single_span = seeds_span.subspan(seed_index, 1)
            | std::views::transform(&search::seed::sequence);

        std::vector<anchor_group> anchor_groups{};
        size_t total_num_raw_anchors = 0;

        // search_n here searches at most max_num_anchors_hard many anchors
        fmindex_collection::search_ng21::search_n(
            index,
            seed_single_span,
            search_scheme,
            config.anchor_choice_strategy == anchor_choice_strategy_t::first_reported ?
                config.max_num_anchors_soft
                : std::max(config.max_num_anchors_hard, config.max_num_anchors_hard + 1),
            [&anchor_groups, &total_num_raw_anchors] (
                [[maybe_unused]] size_t const _seed_index_in_wrapper_range,
                auto cursor,
                size_t const errors
            ) {
                total_num_raw_anchors += cursor.count();
                anchor_groups.emplace_back(cursor, errors);
            }
        );

        if (
            total_num_raw_anchors > config.max_num_anchors_hard
            && config.anchor_choice_strategy != anchor_choice_strategy_t::first_reported
        ) {
            anchors_by_seed.emplace_back(search_result::anchors_of_seed{
                .num_kept_useful_anchors = 0,
                .num_kept_raw_anchors = 0,
                .num_excluded_raw_anchors_by_soft_cap = 0,
                .anchors_by_reference{}
            });

            continue;
        }

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
            case anchor_group_order_t::none:
                break;

            default:
                throw std::runtime_error("(Should be unreachable) internal bug in anchor group order config.");
        }

        // At some point in the future, if the locate operation start showing up in profiles,
        // it could be an optimization opportunity to merge the fmindex cursors before locating

        size_t num_kept_raw_anchors = 0;
        std::vector<anchors_t> anchors_by_reference(num_reference_sequences);
        size_t anchor_group_index = 0;

        // switch case didn't work here, not sure why
        if (config.anchor_choice_strategy == anchor_choice_strategy_t::round_robin) {
            // this is a somewhat complicated implementation using std::set to make sure that
            // the running time is not quadratic in the number of anchor groups
            std::ranges::iota_view init{0ul, anchor_groups.size()};
            std::set<size_t> remaining_group_indices(init.begin(), init.end());
            auto remaining_group_indices_iter = remaining_group_indices.begin();
            size_t round = 0;

            while (
                num_kept_raw_anchors != config.max_num_anchors_soft &&
                !remaining_group_indices.empty()
            ) {
                auto const& [cursor, num_errors] = anchor_groups[*remaining_group_indices_iter];
                // this assumes that cursors are not empty in the beginning
                auto const [reference_id, position] = index.locate(cursor.lb + round);
                anchors_by_reference[reference_id].emplace_back(anchor_t {
                    .pex_leaf_index = seed.pex_leaf_index,
                    .reference_id = reference_id,
                    .reference_position = position,
                    .num_errors = num_errors
                });
                ++num_kept_raw_anchors;

                auto previous_iter = remaining_group_indices_iter;
                ++remaining_group_indices_iter;
                if (cursor.len == round + 1) {
                    remaining_group_indices.erase(previous_iter);
                }

                if (remaining_group_indices_iter == remaining_group_indices.end()) {
                    remaining_group_indices_iter = remaining_group_indices.begin();
                    ++round;
                }
            }
        } else if (
            config.anchor_choice_strategy == anchor_choice_strategy_t::full_groups
            || config.anchor_choice_strategy == anchor_choice_strategy_t::first_reported
        ) {
            while (
                num_kept_raw_anchors != config.max_num_anchors_soft &&
                anchor_group_index < anchor_groups.size()
            ) {
                auto const& [cursor, num_errors] = anchor_groups[anchor_group_index];

                for (auto const& anchor: cursor) {
                    auto const [reference_id, position] = index.locate(anchor);
                    anchors_by_reference[reference_id].emplace_back(anchor_t {
                        .pex_leaf_index = seed.pex_leaf_index,
                        .reference_id = reference_id,
                        .reference_position = position,
                        .num_errors = num_errors
                    });

                    ++num_kept_raw_anchors;
                    if (num_kept_raw_anchors == config.max_num_anchors_soft) {
                        break;
                    }
                }

                ++anchor_group_index;
            }
        } else {
            throw std::runtime_error("(Should be unreachable) internal bug in anchor choice strategy config.");
        }

        size_t const num_excluded_raw_anchors_by_soft_cap = total_num_raw_anchors - num_kept_raw_anchors;

        size_t num_kept_useful_anchors = num_kept_raw_anchors;

        if (config.erase_useless_anchors) {
            num_kept_useful_anchors = erase_useless_anchors(anchors_by_reference);
        }

        anchors_by_seed.emplace_back(search_result::anchors_of_seed{
            .num_kept_useful_anchors = num_kept_useful_anchors,
            .num_kept_raw_anchors = num_kept_raw_anchors,
            .num_excluded_raw_anchors_by_soft_cap = num_excluded_raw_anchors_by_soft_cap,
            .anchors_by_reference = std::move(anchors_by_reference)
        });
    }

    return search_result {
        .anchors_by_seed = std::move(anchors_by_seed),
        .num_fully_excluded_seeds = num_fully_excluded_seeds
    };
}

namespace internal {

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

size_t erase_useless_anchors(std::vector<anchors_t>& anchors_by_reference) {
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

} // namespace internal

} // namespace search
