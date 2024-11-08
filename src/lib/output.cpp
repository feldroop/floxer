#include <math.hpp>
#include <output.hpp>

#include <fstream>
#include <iostream>
#include <limits>
#include <ranges>
#include <stdexcept>

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <ivsigma/ivsigma.h>

#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/std.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/spdlog.h>

namespace output {

void save_index(fmindex const& index, std::filesystem::path const& index_path) {
    spdlog::info("saving index to {}", index_path);

    try {
        auto ofs     = std::ofstream(index_path, std::ios::binary);
        auto archive = cereal::BinaryOutputArchive{ofs};
        archive(index);
    } catch (std::exception const& e) {
        spdlog::warn(
            "An error occured while trying to write the index to "
            "the file {}.\nContinuing without saving the index.\n{}\n",
            index_path,
            e.what()
        );
    }
}

alignment_output::alignment_output(
        std::filesystem::path const& output_path,
        std::vector<input::reference_record> const& references_
) : out(internal::create_seqan_alignment_output(output_path, references_)), references(references_)
{}

void alignment_output::write_alignments_for_query(
    input::query_record const& query,
    alignment::query_alignments alignments
) {
    static constexpr uint8_t mapq_not_available = 255;

    bool primary_alignment_was_written = false;

    for (size_t reference_id = 0; reference_id < references.size(); ++reference_id) {
        auto& reference_alignments = alignments.to_reference(reference_id);
        auto const& reference = references[reference_id];

        for (auto& alignment : reference_alignments) {
            auto flag = alignment.orientation == alignment::query_orientation::reverse_complement ?
                seqan3::sam_flag::on_reverse_strand :
                seqan3::sam_flag::none;

            bool const is_primary_alignment = !primary_alignment_was_written &&
                alignments.best_num_errors().value() == alignment.num_errors;

            std::string query_char_sequence{};
            if (is_primary_alignment) {
                query_char_sequence = ivs::convert_rank_to_char<ivs::d_dna5>(query.rank_sequence);
                primary_alignment_was_written = true;
            } else {
                flag |= seqan3::sam_flag::secondary_alignment;
            }

            seqan3::sam_tag_dictionary tags{};
            using namespace seqan3::literals;
            tags.get<"NM"_tag>() = alignment.num_errors;

            out.emplace_back(
                query.id, // id
                flag, // flag
                reference.id, // ref_id
                math::saturate_value_to_int32_max(alignment.start_in_reference), // ref_offset
                mapq_not_available, // mapq
                std::move(alignment.cigar), // cigar
                query_char_sequence, // seq
                is_primary_alignment ? query.quality : std::string{}, // qual
                tags // tags contains edit distance tag (NM)
            );
        }
    }

    if (!primary_alignment_was_written) {
        out.emplace_back(
            query.id, // id
            seqan3::sam_flag::unmapped, // flag
            std::string{}, // ref_id
            0, // ref_offset
            mapq_not_available, // mapq
            std::vector<seqan3::cigar>{}, // cigar
            ivs::convert_rank_to_char<ivs::d_dna5>(query.rank_sequence), // seq
            query.quality, // qual
            seqan3::sam_tag_dictionary{} // tags
        );
    }
}

void initialize_logger(std::optional<std::filesystem::path> const logfile_path) {
    try {
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
    } catch (std::exception const& e) {
        fmt::print(
            "[ERROR] An error occured while trying to set up logging. "
            "Trying to continue without logging.\n {}\n",
            e.what()
        );
    }
}

std::string format_elapsed_time(std::chrono::duration<double> const elapsed_seconds) {
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

std::string format_large_number(size_t const number) {
    static const char separator = ',';
    static const size_t block_size = 3;

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

    return formatted_number_string;
}

namespace internal {

seqan_alignment_output create_seqan_alignment_output(
    std::filesystem::path const& output_path,
    std::vector<input::reference_record> const& references
) {
    auto reference_id_view = references | std::views::transform(&input::reference_record::id);
    std::vector<std::string> reference_ids(reference_id_view.begin(), reference_id_view.end());

    return seqan3::sam_file_output(
        output_path,
        std::move(reference_ids),
        references | std::views::transform([] (input::reference_record const& ref) {
            return ref.rank_sequence.size();
        }),
        alignment_output_fields_t{}
    );
}

} // namespace internal

} // namespace output
