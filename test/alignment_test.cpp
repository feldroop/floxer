#include <alignment.hpp>

#include <seqan3/io/sam_file/detail/cigar.hpp>

#include <gtest/gtest.h>

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
    std::string const text("ACGACTACTACGAAATTTAAGTATAGGCTACTTTCCGTACGTACGTACGT");
    std::string const pattern("AATTTAAGTCTAGGCTACTTTC");

    auto const text_ptr = reinterpret_cast<const char*>(text.data());
    auto const text_len = static_cast<int>(text.size());

    auto const pattern_ptr = reinterpret_cast<const char*>(pattern.data());
    auto const pattern_len = static_cast<int>(pattern.size());

    int const text_surplus_size = text_len - pattern_len;

    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;

    attributes.distance_metric = distance_metric_t::edit;
    attributes.memory_mode = wavefront_memory_t::wavefront_memory_low;
    attributes.heuristic.strategy = wf_heuristic_strategy::wf_heuristic_none;

    attributes.alignment_form.span = alignment_span_t::alignment_endsfree;
    attributes.alignment_form.pattern_begin_free = 0;
    attributes.alignment_form.pattern_end_free = 0;

    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);

    wf_aligner->alignment_form.text_begin_free = text_surplus_size;
    wf_aligner->alignment_form.text_end_free = text_surplus_size;

    wavefront_align(
        wf_aligner,
        pattern_ptr,
        pattern_len,
        text_ptr,
        text_len
    );

    EXPECT_EQ(wf_aligner->align_status.status, 0);
    EXPECT_EQ(wf_aligner->cigar->score, 1);

    int const alignment_length = wf_aligner->cigar->end_offset -
        wf_aligner->cigar->begin_offset;

    EXPECT_EQ(alignment_length, 50);
    // EXPECT_EQ(alignment_length, 22);

    std::string cigar_buffer(44, '_');
    bool const show_mismatches = true;
    const int written_size = cigar_sprint_SAM_CIGAR(
        cigar_buffer.data(),
        wf_aligner->cigar,
        show_mismatches
    );
    cigar_buffer.resize(written_size);

    std::string expected_cigar("13I9=1X12=15I");
    //std::string expected_cigar("9=1X12=");

    EXPECT_EQ(cigar_buffer, expected_cigar);

    wavefront_aligner_delete(wf_aligner);
}
