#include <fmindex.hpp>
#include <search.hpp>

#include <algorithm>
#include <limits>
#include <ranges>

#include <fmindex-collection/search/SearchNg21.h>
#include <search_schemes/generator/optimum.h>
#include <search_schemes/expand.h>

namespace search {

search_schemes::Scheme const& search_scheme_cache::get(
    size_t const pex_leaf_query_length,
    size_t const pex_leaf_num_errors
) {
    auto const scheme_data = std::make_tuple(pex_leaf_query_length, pex_leaf_num_errors);
    auto const iter = schemes.find(scheme_data);

    if (iter == schemes.end()) {
        auto search_scheme = search_schemes::expand(
            search_schemes::generator::optimum(0, pex_leaf_num_errors),
            pex_leaf_query_length
        );

        auto const [emplaced_iter, _] = schemes.emplace(scheme_data, std::move(search_scheme));
        return emplaced_iter->second;
    }

    return iter->second;
}

bool anchor_t::is_better_than(anchor_t const& other) {
    size_t const position_difference = position < other.position ?
        other.position - position : position - other.position;

    return num_errors <= other.num_errors &&
        position_difference <= other.num_errors - num_errors;
}

void anchor_t::mark_for_erasure() {
    num_errors = internal::erase_marker;
}

bool anchor_t::should_be_erased() const {
    return num_errors == internal::erase_marker;
}

search_result searcher::search_seeds(
    std::vector<seed> const& seeds
) const {
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

        cursors_t cursors(config);

        fmindex_collection::search_ng21::search(
            index,
            seed_single_span,
            search_scheme,
            [&cursors]
            ([[maybe_unused]] size_t const seed_id, auto cursor, size_t const errors) {
                cursors.total_num_raw_anchors += cursor.count();

                auto& cursors_with_given_num_errors = cursors.cursors_by_num_errors.at(errors);
                cursors_with_given_num_errors.total_num_raw_anchors += cursor.count();
                cursors_with_given_num_errors.cursors.emplace_back(cursor);
            }
        );

        size_t num_errors_threshold = 0; // exclude anchors with this number of errors and higher
        size_t num_raw_anchors_below_threshold = 0;

        while (num_errors_threshold <= config.max_num_errors) {
            num_raw_anchors_below_threshold += cursors.cursors_by_num_errors.at(num_errors_threshold).total_num_raw_anchors;

            if (num_raw_anchors_below_threshold >= config.max_num_raw_anchors) {
                break;
            }

            ++num_errors_threshold;
        }

        if (num_errors_threshold == 0) {
            anchors_by_seed.emplace_back(search_result::anchors_of_seed {
                .status = seed_status::fully_excluded,
                .num_kept_useful_anchors = 0,
                .num_excluded_raw_anchors = cursors.total_num_raw_anchors,
                .anchors_by_reference = std::vector<anchors>{}
            });

            ++num_fully_excluded_seeds;
            continue;
        }

        size_t num_kept_raw_anchors = 0;

        std::vector<anchors> anchors_by_reference(num_reference_sequences);
        for (size_t num_errors = 0; num_errors < num_errors_threshold; ++num_errors) {
            auto& cursors_with_given_num_errors = cursors.cursors_by_num_errors.at(num_errors);
            num_kept_raw_anchors += cursors_with_given_num_errors.total_num_raw_anchors;

            for (auto& cursor : cursors_with_given_num_errors.cursors) {
                for (auto const& anchor: cursor) {
                    auto const [reference_id, position] = index.locate(anchor);
                    anchors_by_reference[reference_id].emplace_back(anchor_t {
                        .position = position,
                        .num_errors = num_errors
                    });
                }
            }
        }

        size_t num_useful_anchors = 0;
        for (auto& anchors_of_seed_and_reference : anchors_by_reference) {
            internal::erase_useless_anchors(anchors_of_seed_and_reference);
            num_useful_anchors += anchors_of_seed_and_reference.size();
        }

        size_t const num_excluded_raw_anchors = cursors.total_num_raw_anchors - num_kept_raw_anchors;

        anchors_by_seed.emplace_back(search_result::anchors_of_seed{
            .status = num_excluded_raw_anchors == 0 ? seed_status::not_excluded : seed_status::partly_excluded,
            .num_kept_useful_anchors = num_useful_anchors,
            .num_excluded_raw_anchors = num_excluded_raw_anchors,
            .anchors_by_reference = std::move(anchors_by_reference)
        });
    }

    return search_result {
        .anchors_by_seed = std::move(anchors_by_seed),
        .num_fully_excluded_seeds = num_fully_excluded_seeds
    };
}

cursors_t::cursors_t(search_config const& config) :
    cursors_by_num_errors(config.max_num_errors) {}

namespace internal {

void erase_useless_anchors(anchors& anchors_of_seed_and_reference) {
    // this must stay, otherwise the expression in the below for loop head could underflow
    if (anchors_of_seed_and_reference.empty()) {
        return;
    }

    std::ranges::sort(anchors_of_seed_and_reference, {}, [] (anchor_t const& h) { return h.position; });

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
}

} // namespace internal

} // namespace search
