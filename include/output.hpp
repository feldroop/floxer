#pragma once

#include <about_floxer.hpp>
#include <alignment.hpp>
#include <fmindex.hpp>
#include <input.hpp>

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <ostream>

namespace output {

void save_index(fmindex const& _index, std::filesystem::path const& _index_path);

class sam_output {
private:
    std::ofstream output_stream;

public:
    sam_output(
        std::filesystem::path const& output_path,
        std::vector<input::reference_record> const& reference_records,
        std::string const command_line_call
    );

    void output_for_query(
        input::query_record const& fastq_query,
        std::vector<input::reference_record> const& references,
        alignment::query_alignments const& alignments
    );
};

void initialize_logger(std::optional<std::filesystem::path> const logfile_path);

std::string format_elapsed_time(std::chrono::duration<double> const elapsed_seconds);

std::string format_large_numer(size_t const number, std::string const& unit);

// not thread safe, should only be used once and without other output in between
struct progress_bar{
    size_t const num_updates = 100;
    size_t const max_bar_width = 120;

    char const range_open = '[';
    char const range_close = ']';
    char const bar_char = '=';
    char const bar_tip = '>';
    char const empty_char = ' ';

    size_t const total_num_events;
    size_t next_print_event_index = 0;

    void start();

    void progress(size_t const event_index);

    void finish();

private:
    void print_bar(
        size_t const done_bar_width,
        size_t const remaining_bar_width,
        size_t const percent_done
    );
};

namespace internal {

struct sam_header {
    struct file_level_metadata {
        std::string const version = "1.6"; // VN
        // no alignment sorting order, because floxer only does grouping
        std::string const alignment_grouping = "query"; // GO
    };

    struct reference_sequence_metadata {
        std::string const name; // SN
        uint32_t const length; // LN
    };

    struct program {
        size_t const id = 0; // ID
        std::string const name = about_floxer::program_name; // PN
        std::string const command_line_call; // CL
        std::string const description = about_floxer::short_description + " " + about_floxer::url; // DS
        std::string const version = about_floxer::version; // VN

        program(std::string const command_line_call_);
    };

    file_level_metadata const general_info{}; // @HD
    std::vector<reference_sequence_metadata> reference_sequences{}; // @SQ
    // no @RG read groups for now
    program const floxer_info; // @PG

    sam_header(
        std::vector<input::reference_record> const& reference_records,
        std::string const command_line_call_
    );
};

template<class Traits>
std::basic_ostream<char, Traits>& operator<<(
    std::basic_ostream<char, Traits>& os,
    sam_header::file_level_metadata const& metadata
) {
    os << "@HD\tVN:" << metadata.version
        << "\tGO:" << metadata.alignment_grouping << '\n';

    return os;
}

template<class Traits>
std::basic_ostream<char, Traits>& operator<<(
    std::basic_ostream<char, Traits>& os,
    sam_header::reference_sequence_metadata const& reference_sequence
) {
    return os << "@SQ\tSN:" << reference_sequence.name
        << "\tLN:" << reference_sequence.length << '\n';
}

template<class Traits>
std::basic_ostream<char, Traits>& operator<<(
    std::basic_ostream<char, Traits>& os,
    sam_header::program const& prog
) {
    return os << "@PG\tID:" << prog.id
        << "\tPN:" << prog.name
        << "\tCL:" << prog.command_line_call
        << "\tDS:" << prog.description
        << "\tVN:" << prog.version << '\n';
}

template<class Traits>
std::basic_ostream<char, Traits>& operator<<(
    std::basic_ostream<char, Traits>& os,
    sam_header const& header
) {
    os << header.general_info;
    for (auto const& reference_sequence : header.reference_sequences) {
        os << reference_sequence;
    }
    return os << header.floxer_info;
}

static inline constexpr uint8_t mapq_not_available_marker = 255u;
static inline constexpr std::string string_field_not_available_marker = "*";
static inline constexpr int32_t int_field_not_available_marker = 0;
static inline constexpr int64_t edit_distance_not_available_marker = -1;

struct sam_alignment {
    struct info_flag {
        uint32_t raw_value = 0u;
    };

    // ----- flags that are not used by floxer are commented out -----
    // static inline constexpr info_flag multiple_segments = info_flag { .raw_value = 1u };
    static inline constexpr info_flag each_segment_properly_aligned = info_flag { .raw_value = 2u };
    static inline constexpr info_flag unmapped = info_flag { .raw_value = 4u };
    // static inline constexpr info_flag next_unmapped = info_flag { .raw_value = 8u };
    static inline constexpr info_flag seq_reverse_complemented = info_flag { .raw_value = 16u };
    // static inline constexpr info_flag next_seq_reverse_complemented = info_flag { .raw_value = 32u };
    static inline constexpr info_flag first_segment = info_flag { .raw_value = 64u };
    static inline constexpr info_flag last_segment = info_flag { .raw_value = 128u };
    static inline constexpr info_flag secondary_alignment = info_flag { .raw_value = 256u };
    // static inline constexpr info_flag not_passing_filters = info_flag { .raw_value = 512u };
    // static inline constexpr info_flag pcr_or_optical_duplicate = info_flag { .raw_value = 1024u };
    // static inline constexpr info_flag supplementary_alignment = info_flag { .raw_value = 2048u };

    std::string const& qname;
    info_flag const flag;
    std::string const& rname;
    int32_t const pos;
    uint16_t const mapq;
    std::string const& cigar;
    std::string const& rnext;
    int32_t const pnext;
    int32_t const tlen;
    std::string const& seq;
    std::string const& qual;
    int64_t const custom_field_edit_distance;
};

sam_alignment::info_flag operator|(
    sam_alignment::info_flag const& left,
    sam_alignment::info_flag const& right
);

sam_alignment::info_flag& operator|=(
    sam_alignment::info_flag& left,
    sam_alignment::info_flag const& right
);

template<class Traits>
std::basic_ostream<char, Traits>& operator<<(
    std::basic_ostream<char, Traits>& os,
    sam_alignment const& alignment
) {
    return os << alignment.qname << '\t'
        << alignment.flag.raw_value << '\t'
        << alignment.rname << '\t'
        << alignment.pos << '\t'
        << alignment.mapq << '\t'
        << alignment.cigar << '\t'
        << alignment.rnext << '\t'
        << alignment.pnext << '\t'
        << alignment.tlen << '\t'
        << alignment.seq << '\t'
        << alignment.qual << '\t'
        << "NM:i:" << alignment.custom_field_edit_distance << '\n';
}

int32_t saturate_value_to_int32_max(size_t const value);

} // namespace internal

} // namespace output
