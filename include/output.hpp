#pragma once

#include <about_floxer.hpp>
#include <alignment.hpp>
#include <fmindex.hpp>
#include <input.hpp>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>

#include <spdlog/stopwatch.h>

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

std::string format_elapsed_time(spdlog::stopwatch const& stopwatch);

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

} // namespace output
