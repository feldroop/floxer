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

using alignment_output_t = seqan3::sam_file_output<
    alignment_output_fields_t,
    seqan3::type_list<seqan3::format_sam, seqan3::format_bam>,
    std::vector<std::string>
>;

alignment_output_t create_alignment_output(
    std::filesystem::path const& output_path,
    std::vector<input::reference_record> const& references
);

void output_for_query(
    alignment_output_t& alignment_output,
    input::query_record const& fastq_query,
    std::vector<input::reference_record> const& references,
    alignment::query_alignments alignments
);

void initialize_logger(std::optional<std::filesystem::path> const logfile_path);

std::string format_elapsed_time(std::chrono::duration<double> const elapsed_seconds);

std::string format_large_numer(size_t const number);

} // namespace output
