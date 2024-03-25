#include <output.hpp>

#include <fstream>
#include <limits>
#include <ostream>

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <fmt/core.h>

namespace output {

void save_index(fmindex const& _index, std::filesystem::path const& _index_path) {
    auto ofs     = std::ofstream(_index_path, std::ios::binary);
    auto archive = cereal::BinaryOutputArchive{ofs};
    archive(_index);
}

int32_t saturate_value_to_int32_max(size_t const value) {
    if (value > std::numeric_limits<int32_t>::max()) {
        return std::numeric_limits<uint32_t>::max();
    } else {
        return static_cast<uint32_t>(value);
    }
}

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
        // no command line call, because it might leak unexpected things
        std::string const description = about_floxer::long_description + " " + about_floxer::url; // DS
        std::string const version = about_floxer::version; // VN
    };

    file_level_metadata const general_info{}; // @HD
    std::vector<reference_sequence_metadata> reference_sequences{}; // @SQ
    // no @RG read groups for now
    program const floxer_info{}; // @PG
    std::string const comment_line; // @CO

    sam_header(
        std::vector<input::reference_record> const& reference_records,
        std::string const comment_line_
    ) : comment_line{std::move(comment_line_)} {
        for (auto const& record : reference_records) {
            uint32_t reference_length = saturate_value_to_int32_max(record.sequence_length);
            if (record.sequence_length > std::numeric_limits<int32_t>::max()) {
                fmt::println(
                    "[OUTPUT WARNING]: the sequence {} is too long for the SAM file format (length {})\n"
                    "Values in the output file will be set to INT32_MAX.",
                    record.raw_tag,
                    record.sequence_length
                );
            }
            
            reference_sequences.push_back(reference_sequence_metadata {
                .name = record.sam_format_sanitized_name,
                .length = reference_length
            });
        }
    }
};

template<class Traits>
std::basic_ostream<char, Traits>& operator<<(
    std::basic_ostream<char, Traits>& os,
    sam_header::file_level_metadata const& metadata
) {
    os << "@HG\tVN:" << metadata.version
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
    os << header.floxer_info;
    return os << "@CO\t" << header.comment_line;
}

static constexpr uint8_t mapq_not_available_marker = 255u;
static constexpr std::string string_field_not_available_marker = "*";
static constexpr int32_t int_field_not_available_marker = 0;
static constexpr uint64_t edit_distance_not_available_marker = 999999ul;

struct sam_alignment {
    struct info_flag {
        uint32_t raw_value = 0u;
    };

    // ----- flags that are not used by floxer are commented out -----
    // static constexpr info_flag multiple_segments = info_flag { .raw_value = 1u };
    static constexpr info_flag each_segment_properly_aligned = info_flag { .raw_value = 2u };
    static constexpr info_flag unmapped = info_flag { .raw_value = 4u };
    // static constexpr info_flag next_unmapped = info_flag { .raw_value = 8u };
    static constexpr info_flag seq_reverse_complemented = info_flag { .raw_value = 16u };
    // static constexpr info_flag next_seq_reverse_complemented = info_flag { .raw_value = 32u };
    static constexpr info_flag first_segment = info_flag { .raw_value = 64u };
    static constexpr info_flag last_segment = info_flag { .raw_value = 128u };
    static constexpr info_flag secondary_alignment = info_flag { .raw_value = 256u };
    // static constexpr info_flag not_passing_filters = info_flag { .raw_value = 512u };
    // static constexpr info_flag pcr_or_optical_duplicate = info_flag { .raw_value = 1024u };
    // static constexpr info_flag supplementary_alignment = info_flag { .raw_value = 2048u };

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
    uint64_t const custom_field_edit_distance;
};

sam_alignment::info_flag operator|(
    sam_alignment::info_flag const& left,
    sam_alignment::info_flag const& right
) {
    return sam_alignment::info_flag{ .raw_value = left.raw_value | right.raw_value };
}

sam_alignment::info_flag& operator|=(
    sam_alignment::info_flag& left,
    sam_alignment::info_flag const& right
) {
    left.raw_value |= right.raw_value;
    return left;
}

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

sam_output::sam_output(
    std::filesystem::path const& output_path,
    std::vector<input::reference_record> const& reference_records,
    std::string const comment_line
) : output_stream(output_path) {
    output_stream << sam_header(reference_records, std::move(comment_line));
}

void sam_output::output_for_query(
    input::query_record const& fastq_query,
    std::vector<input::reference_record> const& references,
    verification::fastq_query_alignments const& alignments
) {
    bool found_any_alignments = false;

    for (size_t reference_id = 0; reference_id < references.size(); ++reference_id) {
        auto const& reference_alignments = alignments.for_reference(reference_id);
        auto const& reference = references[reference_id];
    
        if (!reference_alignments.empty()) {
            found_any_alignments = true;
        }

        for (auto const& alignment : std::views::values(reference_alignments)) {
            sam_alignment::info_flag flag = sam_alignment::each_segment_properly_aligned
                | sam_alignment::first_segment
                | sam_alignment::last_segment;

            bool const is_primary_alignment = alignments.is_primary_alignment(alignment);
            if (!is_primary_alignment) {
                flag |= sam_alignment::secondary_alignment;
            }

            if (alignment.is_reverse_complement) {
                flag |= sam_alignment::seq_reverse_complemented;
            }

            output_stream << sam_alignment {
                .qname = fastq_query.sam_format_sanitized_name, // floxer assumes single segment template
                .flag = flag.raw_value,
                .rname = reference.sam_format_sanitized_name,
                .pos = saturate_value_to_int32_max(alignment.start_in_reference + 1), // .sam format is 1-based here
                .mapq = mapq_not_available_marker,
                .cigar = alignment.cigar.to_string(),
                .rnext = string_field_not_available_marker, // floxer assumes single segment template
                .pnext = int_field_not_available_marker, // floxer assumes single segment template
                .tlen = int_field_not_available_marker, // floxer assumes single segment template
                .seq = is_primary_alignment ? fastq_query.char_sequence : string_field_not_available_marker,
                .qual = (is_primary_alignment && !fastq_query.quality.empty()) ? 
                    fastq_query.quality : string_field_not_available_marker,
                .custom_field_edit_distance = alignment.num_errors
            };
        }
    }

    if (!found_any_alignments) {
        sam_alignment::info_flag flag = sam_alignment::unmapped
            | sam_alignment::first_segment
            | sam_alignment::last_segment;

        output_stream << sam_alignment {
            .qname = fastq_query.sam_format_sanitized_name, // floxer assumes single segment template
            .flag = flag.raw_value,
            .rname = string_field_not_available_marker,
            .pos = int_field_not_available_marker,
            .mapq = mapq_not_available_marker,
            .cigar = string_field_not_available_marker,
            .rnext = string_field_not_available_marker, // floxer assumes single segment template
            .pnext = int_field_not_available_marker, // floxer assumes single segment template
            .tlen = int_field_not_available_marker, // floxer assumes single segment template
            .seq = fastq_query.char_sequence,
            .qual = fastq_query.quality.empty() ? string_field_not_available_marker: fastq_query.quality,
            .custom_field_edit_distance = edit_distance_not_available_marker
        };
    }
}

} // namespace output
