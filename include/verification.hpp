#pragma once

#include <alignment.hpp>
#include <input.hpp>
#include <intervals.hpp>
#include <pex.hpp>
#include <search.hpp>
#include <statistics.hpp>

namespace verification {

namespace internal {

struct span_config;

}

struct query_verifier {
    // should only be called once on each instance
    void verify(pex::verification_kind_t const kind);

private:
    void direct_full_verification();

    void hierarchical_verification();

    bool root_was_already_verified() const;

    internal::span_config compute_root_reference_span_config() const;

public:
    pex::pex_tree const& pex_tree;
    search::anchor_t const& anchor;
    pex::pex_tree::node pex_node;
    std::span<const uint8_t> const query;
    alignment::query_orientation const orientation;
    input::reference_record const& reference;
    intervals::verified_intervals& already_verified_intervals;
    double const extra_verification_ratio;
    alignment::query_alignments& alignments;
    statistics::search_and_alignment_statistics& stats;
};

namespace internal {

struct span_config {
    size_t const offset{};
    size_t const length{};

    size_t const applied_extra_verification_length_per_side;

    intervals::half_open_interval as_half_open_interval() const;
};

span_config compute_reference_span_start_and_length(
    search::anchor_t const& anchor,
    pex::pex_tree::node const& pex_node,
    size_t const leaf_query_index_from,
    size_t const full_reference_length,
    double const extra_verification_ratio
);

alignment::alignment_outcome try_to_align_pex_node_query_with_reference_span(
    pex::pex_tree::node const& pex_node,
    input::reference_record const& reference,
    span_config const reference_span_config,
    std::span<const uint8_t> const query,
    alignment::query_orientation const orientation,
    alignment::query_alignments& alignments,
    statistics::search_and_alignment_statistics& stats
);

}

} // verification
