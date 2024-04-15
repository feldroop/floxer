#include <output.hpp>

#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <ranges>
#include <stdexcept>

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/std.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/spdlog.h>

#include <ivsigma/ivsigma.h>

namespace output {

void save_index(fmindex const& _index, std::filesystem::path const& _index_path) {
    auto ofs     = std::ofstream(_index_path, std::ios::binary);
    auto archive = cereal::BinaryOutputArchive{ofs};
    archive(_index);
}

sam_output::sam_output(
    std::filesystem::path const& output_path,
    std::vector<input::reference_record> const& reference_records,
    std::string const command_line_call
) : output_stream(output_path) {
    output_stream << internal::sam_header(reference_records, std::move(command_line_call));
}

void sam_output::output_for_query(
    input::query_record const& fastq_query,
    std::vector<input::reference_record> const& references,
    alignment::query_alignments const& alignments
) {
    using namespace internal;

    bool found_any_alignments = false;

    for (size_t reference_id = 0; reference_id < references.size(); ++reference_id) {
        auto const& reference_alignments = alignments.to_reference(reference_id);
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
                .seq = is_primary_alignment ? ivs::convert_rank_to_char<ivs::d_dna5>(fastq_query.rank_sequence)
                    : string_field_not_available_marker,
                .qual = (is_primary_alignment && !fastq_query.quality.empty()) ? 
                    fastq_query.quality : string_field_not_available_marker,
                .custom_field_edit_distance = static_cast<int64_t>(alignment.num_errors)
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
            .seq = ivs::convert_rank_to_char<ivs::d_dna5>(fastq_query.rank_sequence),
            .qual = fastq_query.quality.empty() ? string_field_not_available_marker: fastq_query.quality,
            .custom_field_edit_distance = edit_distance_not_available_marker
        };
    }
}

void initialize_logger(std::optional<std::filesystem::path> const logfile_path) {
    std::vector<spdlog::sink_ptr> sinks;
    
    auto console_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
    console_sink->set_level(spdlog::level::info);

    sinks.push_back(console_sink);

    if (logfile_path.has_value()) {
        auto const max_logfile_size = 1024 * 1024 * 5; // 5 MB
        auto const max_num_logfiles = 1;

        auto file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
            logfile_path.value(),
            max_logfile_size,
            max_num_logfiles
        );
        file_sink->set_level(spdlog::level::trace);
        std::string const debug_log_pattern = "[thread %t] %+";
        file_sink->set_pattern(debug_log_pattern);

        sinks.push_back(file_sink);
    }

    auto logger = std::make_shared<spdlog::logger>(about_floxer::program_name, begin(sinks), end(sinks));
    logger->set_level(spdlog::level::trace);
    logger->flush_on(spdlog::level::debug);
    
    spdlog::set_default_logger(logger);
}

std::string format_elapsed_time(spdlog::stopwatch const& stopwatch) {
    auto const elapsed_seconds = stopwatch.elapsed();
    if (elapsed_seconds <= std::chrono::seconds(60)) {
        return fmt::format("{:.7} seconds", elapsed_seconds);
    }

    size_t const all_in_seconds = elapsed_seconds.count();
    size_t const seconds = all_in_seconds % 60;

    size_t const all_in_minutes = all_in_seconds / 60;
    size_t const minutes = all_in_minutes % 60;

    size_t const all_in_hours = all_in_minutes / 60;
    size_t const hours = all_in_hours % 24;

    if (hours > 0) {
        return fmt::format("{}:{:02}:{:02} hours", hours, minutes, seconds);
    } else {
        return fmt::format("{:02}:{:02} minutes", minutes, seconds);
    }
}

void progress_bar::progress(size_t const event_index) {
    if (event_index >= next_print_event_index) {
        double const fraction_done = event_index / static_cast<double>(total_num_events);
        size_t const done_bar_width = max_bar_width * fraction_done;
        size_t const remaining_bar_width = max_bar_width - done_bar_width;
        size_t const percent_done = fraction_done * 100;

        print_bar(done_bar_width, remaining_bar_width, percent_done);

        next_print_event_index += total_num_events / num_updates;
    }
}

void progress_bar::start() {
    print_bar(0, max_bar_width, 0);
}

void progress_bar::finish() {
    print_bar(max_bar_width, 0, 100);
    std::cerr << '\n';
}

void progress_bar::print_bar(
    size_t const done_bar_width,
    size_t const remaining_bar_width,
    size_t const percent_done
) {
    std::string bar{};
    bar += range_open;
    bar += std::string(done_bar_width, bar_char);
    bar += bar_tip;
    bar += std::string(remaining_bar_width, empty_char);
    bar += range_close;

    fmt::print(stderr, "\rProgress: {} {: >3}%", bar, percent_done);
    std::cerr.flush();
}

std::string format_large_numer(size_t const number, std::string const& unit) {
    static const char separator = ' ';
    static const size_t block_size = 3;
    static const size_t unit_prefix_base = 1000;
    static std::vector<std::string> const unit_prefixes{"", "kilo", "mega", "giga", "tera", "peta"};

    std::string const raw_number_string = fmt::format("{}", number);
    
    std::string formatted_number_string{};
    size_t i = 0;
    for(auto const digit : std::views::reverse(raw_number_string)) {
        if (i > 0 && i % block_size == 0) {
            formatted_number_string += separator;
        }

        formatted_number_string += digit;
        ++i;
    }
    std::ranges::reverse(formatted_number_string);

    double number_with_unit_prefix = number;
    std::string chosen_unit_prefix{};
    for (auto const& unit_prefix : unit_prefixes) {
        if (number_with_unit_prefix >= unit_prefix_base && unit_prefix != unit_prefixes.back()) {
            number_with_unit_prefix /= unit_prefix_base;
            continue;
        }

        chosen_unit_prefix = unit_prefix;
        break;
    }

    if (number >= 1000) {
        formatted_number_string += fmt::format(
            " ({:.2} {}{})",
            number_with_unit_prefix,
            chosen_unit_prefix,
            unit
        );
    }

    return formatted_number_string;
}

namespace internal {

sam_header::program::program(std::string const command_line_call_)
    : command_line_call{std::move(command_line_call_)} {}

sam_header::sam_header(
        std::vector<input::reference_record> const& reference_records,
        std::string const command_line_call_
    ) : floxer_info(std::move(command_line_call_)) {
    for (auto const& record : reference_records) {
        uint32_t reference_length = saturate_value_to_int32_max(record.rank_sequence.size());
        if (record.rank_sequence.size() > std::numeric_limits<int32_t>::max()) {
            spdlog::warn(
                "The sequence {} is too long for the SAM file format (length {})\n"
                "Values in the output file will be set to INT32_MAX.",
                record.raw_tag,
                record.rank_sequence.size()
            );
        }
        
        reference_sequences.push_back(reference_sequence_metadata {
            .name = record.sam_format_sanitized_name,
            .length = reference_length
        });
    }
}

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

int32_t saturate_value_to_int32_max(size_t const value) {
    if (value > std::numeric_limits<int32_t>::max()) {
        return std::numeric_limits<int32_t>::max();
    } else {
        return static_cast<int32_t>(value);
    }
}

} // namespace internal

} // namespace output
