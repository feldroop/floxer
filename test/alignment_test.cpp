#include <alignment.hpp>

#include <seqan3/io/sam_file/detail/cigar.hpp>

#include <gtest/gtest.h>

TEST(alignment, small_wrapped_seqan3) {
    using namespace alignment;
    set_alignment_backend_global(alignment_backend::seqan3);

    std::vector<uint8_t> reference{ 0, 0, 1, 2, 1, 3, 0, 2, 2, 3, 0, 1 };
    std::vector<uint8_t> query{ 1, 2, 1, 3, 1, 2, 2 };

    alignment_config const config{
        .reference_span_offset = 0,
        .num_allowed_errors = 2,
        .orientation = query_orientation::forward,
        .mode = alignment_mode::verify_and_return_alignment
    };

    aligner seqan3_aligner;
    auto const result = seqan3_aligner.align(reference, query, config);

    EXPECT_EQ(result.outcome, alignment_outcome::alignment_exists);

    EXPECT_TRUE(result.alignment.has_value());

    EXPECT_EQ(result.alignment.value().num_errors, 1);
    EXPECT_EQ(result.alignment.value().orientation, query_orientation::forward);
    EXPECT_EQ(result.alignment.value().start_in_reference, 2);
    EXPECT_EQ(result.alignment.value().cigar, seqan3::detail::parse_cigar("4=1X2="));
}

TEST(alignment, small_wrapped_wfa2) {
    using namespace alignment;
    set_alignment_backend_global(alignment_backend::wfa2);

    std::vector<uint8_t> reference{ 0, 0, 1, 2, 1, 3, 0, 2, 2, 3, 0, 1 };
    std::vector<uint8_t> query{ 1, 2, 1, 3, 1, 2, 2 };

    alignment_config const config{
        .reference_span_offset = 0,
        .num_allowed_errors = 2,
        .orientation = query_orientation::forward,
        .mode = alignment_mode::verify_and_return_alignment
    };

    aligner wfa2_aligner;
    auto const result = wfa2_aligner.align(reference, query, config);

    EXPECT_EQ(result.outcome, alignment_outcome::alignment_exists);

    EXPECT_TRUE(result.alignment.has_value());

    EXPECT_EQ(result.alignment.value().num_errors, 1);
    EXPECT_EQ(result.alignment.value().orientation, query_orientation::forward);
    EXPECT_EQ(result.alignment.value().start_in_reference, 2);
    EXPECT_EQ(result.alignment.value().cigar, seqan3::detail::parse_cigar("4=1X2="));
}

TEST(alignment, small_direct_wfa2) {
    // "glocal" example from WFA2-lib README
    std::string const reference("ACGACTACTACGAAATTTAAGTATAGGCTACTTTCCGTACGTACGTACGT");
    std::string const query    (             "AATTTAAGTCTAGGCTACTTTC"               );

    auto const reference_ptr = reinterpret_cast<const char*>(reference.data());
    auto const reference_len = static_cast<int>(reference.size());

    auto const query_ptr = reinterpret_cast<const char*>(query.data());
    auto const query_len = static_cast<int>(query.size());

    int const reference_size_surplus = reference_len - query_len;

    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;

    attributes.distance_metric = distance_metric_t::edit;
    attributes.memory_mode = wavefront_memory_t::wavefront_memory_low;
    attributes.heuristic.strategy = wf_heuristic_strategy::wf_heuristic_none;

    attributes.alignment_form.span = alignment_span_t::alignment_endsfree;
    // WFA2 text is the query in SAM format
    attributes.alignment_form.text_begin_free = 0;
    attributes.alignment_form.text_end_free = 0;

    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);

    // WFA2 pattern is the reference in SAM format
    wf_aligner->alignment_form.pattern_begin_free = reference_size_surplus;
    wf_aligner->alignment_form.pattern_end_free = reference_size_surplus;

    wavefront_align(
        wf_aligner,
        reference_ptr,
        reference_len,
        query_ptr,
        query_len
    );

    EXPECT_EQ(wf_aligner->align_status.status, 0);
    EXPECT_EQ(wf_aligner->cigar->score, 1);

    // this trims the end insertions but not the ones at the front
    // bool const trimmed = cigar_maxtrim_gap_linear(wf_aligner->cigar, &wf_aligner->penalties.linear_penalties);
    // EXPECT_TRUE(trimmed);

    // these seem to be not the offset regarding ends free gaps, but something else
    // EXPECT_EQ(wf_aligner->cigar->begin_offset, 13);
    // EXPECT_EQ(wf_aligner->cigar->end_offset, 35);

    int const alignment_length = wf_aligner->cigar->end_offset -
        wf_aligner->cigar->begin_offset;

    EXPECT_GE(alignment_length, 0);

    std::string cigar_buffer(alignment_length * 2, '\0');
    bool const show_mismatches = true;
    const int written_size = cigar_sprint_SAM_CIGAR(
        cigar_buffer.data(),
        wf_aligner->cigar,
        show_mismatches
    );
    cigar_buffer.resize(written_size);

    std::string expected_cigar("9=1X12=");

    // due to the unexpected API behavior above, we need to do the endsfree trimming by ourselves
    auto const trim_result = alignment::internal::trim_wfa2_cigar(cigar_buffer);

    EXPECT_EQ(trim_result.cigar, expected_cigar);
    EXPECT_EQ(trim_result.num_trimmed_start, 13);

    wavefront_aligner_delete(wf_aligner);
}

TEST(alignment, trim_wfa2_cigar) {
    auto const trim_result_basic = alignment::internal::trim_wfa2_cigar("428D8=5I3X6D5X87=72D");

    EXPECT_EQ(trim_result_basic.num_trimmed_start, 428);
    EXPECT_EQ(trim_result_basic.cigar, "8=5I3X6D5X87=");

    std::string const fine_cigar = "1=2D5I3X6D5X87=";
    auto const trim_result_nothing_done = alignment::internal::trim_wfa2_cigar(fine_cigar);

    EXPECT_EQ(trim_result_nothing_done.num_trimmed_start, 0);
    EXPECT_EQ(trim_result_nothing_done.cigar, fine_cigar);
}
