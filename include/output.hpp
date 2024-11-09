#pragma once

#include <about_floxer.hpp>
#include <alignment.hpp>
#include <fmindex.hpp>
#include <input.hpp>

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>

#include <seqan3/core/debug_stream/tuple.hpp>
#include <seqan3/io/sam_file/output.hpp>

namespace output {

void save_index(fmindex const& _index, std::filesystem::path const& _index_path);

using alignment_output_fields_t = seqan3::fields<
    seqan3::field::id,
    seqan3::field::flag,
    seqan3::field::ref_id,
    seqan3::field::ref_offset,
    seqan3::field::mapq,
    seqan3::field::cigar,
    seqan3::field::seq,
    seqan3::field::qual,
    seqan3::field::tags
>;

using seqan_alignment_output = seqan3::sam_file_output<
    alignment_output_fields_t,
    seqan3::type_list<seqan3::format_sam, seqan3::format_bam>,
    std::vector<std::string>
>;

// simple wrapper of the seqan3 output to write alignments in the way I want
class alignment_output {
public:
    alignment_output(
        std::filesystem::path const& output_path,
        std::vector<input::reference_record> const& references
    );

    void write_alignments_for_query(
        input::query_record const& fastq_query,
        alignment::query_alignments alignments
    );

private:
    seqan_alignment_output out;
    std::vector<input::reference_record> const& references;
};

void initialize_logger(std::optional<std::filesystem::path> const logfile_path, bool const console_debug_logs);

std::string format_elapsed_time(std::chrono::duration<double> const elapsed_seconds);

std::string format_large_number(size_t const number);

namespace internal {

seqan_alignment_output create_seqan_alignment_output(
    std::filesystem::path const& output_path,
    std::vector<input::reference_record> const& references
);

} // namespace internal

} // namespace output
