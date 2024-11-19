#include <alignment.hpp>
#include <input.hpp>
#include <pex.hpp>
#include <statistics.hpp>
#include <verification.hpp>

#include <gtest/gtest.h>
#include <seqan3/io/sam_file/detail/cigar.hpp>
#include <seqan3/core/debug_stream.hpp>

TEST(verification, verify) {
    input::reference_record const reference {
        .id = "",
        .rank_sequence = {
        4,2,3,4,3,4,4,4,3,2,
        4,3,3,2,2,3,4,4,3,3,
        4,3,2,2,1,4,3,3,4,2,
        4,4,4,3,3,2,1,1,1,2,
        3,4,4,3,2,4,4,2,1,4,
        4,3,4,4,4,4,3,3,2,1, // query
        2,3,4,3,2,1,2,3,4,3, // query
        1,4,2,1,4,4,2,2,3,4, // query
        3,3,2,1,4,4,1,1,1,2,
        4,3,2,1,2,2,2,3,3,1
        },
        .internal_id = 0
    };

    std::vector<uint8_t> query {
        4,3,4,4,4,4,3,3,2,1,4, // insertion at end
        2,3,4,3,2,1,2,3,4, // deletion at end
        1,4,2,1,4,4,2,2,3,4
    };


    pex::pex_tree pex_tree(pex::pex_tree_config{
        query.size(),
        5,
        1,
        pex::pex_tree_build_strategy::bottom_up
    });

    search::anchor_t anchor {
        .pex_leaf_index = 0,
        .reference_id = 0,
        .reference_position = 50,
        .num_errors = 0
    };

    auto const& pex_node = pex_tree.get_leaves().at(0);

    shared_mutex_guarded<intervals::verified_intervals> already_verified_intervals;

    double const extra_verification_ratio = 0.1;

    alignment::query_alignments alignments(1);
    statistics::search_and_alignment_statistics stats;

    verification::query_verifier verifier {
        .pex_tree = pex_tree,
        .anchor = anchor,
        .pex_leaf_node = pex_node,
        .query = query,
        .orientation = alignment::query_orientation::reverse_complement,
        .reference = reference,
        .already_verified_intervals = already_verified_intervals,
        .extra_verification_ratio = extra_verification_ratio,
        .alignments = alignments,
        .stats = stats
    };

    verifier.verify(pex::verification_kind_t::hierarchical);

    EXPECT_EQ(alignments.size(), 1);
    auto const& alignment = alignments.to_reference(0).at(0);

    EXPECT_EQ(alignment.cigar, seqan3::detail::parse_cigar("10=1I9=1D10="));
    EXPECT_EQ(alignment.num_errors, 2);
    EXPECT_EQ(alignment.orientation, alignment::query_orientation::reverse_complement);
    EXPECT_EQ(alignment.start_in_reference, 50);

    verifier.verify(pex::verification_kind_t::hierarchical);

    // nothing should change because of already_verified_intervals
    EXPECT_EQ(alignments.size(), 1);

    shared_mutex_guarded<intervals::verified_intervals> deactivated_already_verified_intervals;
    {
        auto && [lock, ivls] = already_verified_intervals.lock_unique();
        ivls.configure(intervals::use_interval_optimization::off);
    }
    verification::query_verifier direct_verifier {
        .pex_tree = pex_tree,
        .anchor = anchor,
        .pex_leaf_node = pex_node,
        .query = query,
        .orientation = alignment::query_orientation::reverse_complement,
        .reference = reference,
        .already_verified_intervals = deactivated_already_verified_intervals,
        .extra_verification_ratio = extra_verification_ratio,
        .alignments = alignments,
        .stats = stats
    };

    direct_verifier.verify(pex::verification_kind_t::direct_full);

    EXPECT_EQ(alignments.size(), 2);
    EXPECT_EQ(alignments.to_reference(0).at(1), alignments.to_reference(0).at(0));

    // add more errors such that no alignment should be added
    query[5] = 1;
    query[6] = 1;
    query[11] = 3;
    query[20] = 2;

    direct_verifier.verify(pex::verification_kind_t::direct_full);

    EXPECT_EQ(alignments.size(), 2);
}


TEST(verification, compute_reference_span_start_and_length) {
    search::anchor_t anchor {
        .pex_leaf_index = 0,
        .reference_id = 0,
        .reference_position = 100'755,
        .num_errors = 25
    };

    pex::pex_tree::node pex_node {
        .parent_id = 0,
        .query_index_from = 500,
        .query_index_to = 999,
        .num_errors = 30
    };

    size_t const leaf_query_index_from = 750;
    size_t const full_reference_length = 1'000'000;
    double const no_extra_verification_rate = 0.0;

    auto const span_config_base = verification::internal::compute_reference_span_start_and_length(
        anchor, pex_node, leaf_query_index_from, full_reference_length, no_extra_verification_rate
    );

    EXPECT_EQ(span_config_base.offset, 100'475);
    EXPECT_EQ(span_config_base.length, 561);
    EXPECT_EQ(span_config_base.applied_extra_verification_length_per_side, 0);

    double const extra_verification_rate = 0.01;
    auto const span_config_with_extra_space = verification::internal::compute_reference_span_start_and_length(
        anchor, pex_node, leaf_query_index_from, full_reference_length, extra_verification_rate
    );

    EXPECT_EQ(span_config_with_extra_space.offset, 100'469);
    EXPECT_EQ(span_config_with_extra_space.length, 573);
    EXPECT_EQ(span_config_with_extra_space.applied_extra_verification_length_per_side, 6); // ceil(561 * 0.01)
}

TEST(verification, try_to_align_pex_node_query_with_reference_span) {
    pex::pex_tree::node pex_node {
        .parent_id = pex::pex_tree::node::null_id, // fake root node
        .query_index_from = 40,
        .query_index_to = 84,
        .num_errors = 5
    };

    input::reference_record const reference {
        .id = "",
        .rank_sequence = {
            2,2,2,2,2,2,2,2,2,2,
            1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,
            2,2,2,2,2,2,2,2,2,2
        },
        .internal_id = 0
    };

    verification::internal::span_config const span_config {
        .offset = 50,
        .length = 50,
        .applied_extra_verification_length_per_side = 0
    };

    std::vector<uint8_t> query {
        1,1,1,3,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,3,
        1,4,1,1,1,2,1,1,1,1,
        1,1,1,3,1,1,1,4,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1
    };

    alignment::query_alignments alignments(1);
    statistics::search_and_alignment_statistics stats;

    auto const outcome_expected_exists = verification::internal::try_to_align_pex_node_query_with_reference_span(
        pex_node,
        reference,
        span_config,
        std::span(query),
        alignment::query_orientation::forward,
        alignments,
        stats
    );

    EXPECT_EQ(outcome_expected_exists, alignment::alignment_outcome::alignment_exists);

    EXPECT_EQ(alignments.size(), 1);

    auto const& alignment = alignments.to_reference(0).at(0);

    EXPECT_EQ(alignment.num_errors, 5);
    EXPECT_EQ(alignment.orientation, alignment::query_orientation::forward);
    EXPECT_EQ(alignment.start_in_reference, 50);

    pex_node.parent_id = 0; // ---------- not root anymore ----------

    auto const outcome_expected_exists_again = verification::internal::try_to_align_pex_node_query_with_reference_span(
        pex_node,
        reference,
        span_config,
        std::span(query),
        alignment::query_orientation::forward,
        alignments,
        stats
    );

    EXPECT_EQ(outcome_expected_exists_again, alignment::alignment_outcome::alignment_exists);

    EXPECT_EQ(alignments.size(), 1);

    query[42] = 2; // ---------- too many errors ----------
    auto const outcome_expected_does_not_exist = verification::internal::try_to_align_pex_node_query_with_reference_span(
        pex_node,
        reference,
        span_config,
        std::span(query),
        alignment::query_orientation::forward,
        alignments,
        stats
    );

    EXPECT_EQ(outcome_expected_does_not_exist, alignment::alignment_outcome::no_adequate_alignment_exists);
    EXPECT_EQ(alignments.size(), 1);
}
