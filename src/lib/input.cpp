#include <input.hpp>

#include <algorithm>
#include <fstream>
#include <numeric>
#include <ranges>
#include <sstream>
#include <unordered_set>

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>

#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/std.h>
#include <spdlog/spdlog.h>

namespace input {

size_t query_record::num_errors_from_user_config(cli::command_line_input const& cli_input) const {
    if (cli_input.query_error_probability().has_value()) {
        double const num_errors_frac = rank_sequence.size() * cli_input.query_error_probability().value();

        // handle floating point inaccuracy
        static constexpr double epsilon = 0.000000001;
        if (std::abs(num_errors_frac - std::round(num_errors_frac)) < epsilon) {
            return static_cast<size_t>(std::round(num_errors_frac) + epsilon);
        } else {
            return static_cast<size_t>(std::ceil(num_errors_frac));
        }
    } else {
        return cli_input.query_num_errors().value();
    }
}

references read_references(std::filesystem::path const& reference_sequence_path) {
    spdlog::info("reading reference sequences from {}", reference_sequence_path);

    std::vector<reference_record> records{};

    size_t internal_id = 0;
    size_t total_length = 0;

    for (auto const record_view : ivio::fasta::reader{{ .input = reference_sequence_path }}) {
        std::string const id = internal::extract_record_id(record_view.id);

        if (record_view.seq.empty()) {
            spdlog::warn(
                "The record {} in the reference file has an empty sequence and will be skipped.\n",
                id
            );

            continue;
        }

        std::vector<uint8_t> const rank_sequence = internal::chars_to_rank_sequence(record_view.seq);

        spdlog::debug("read reference, id: {}, length {}", id, rank_sequence.size());

        total_length += rank_sequence.size();

        records.emplace_back(
            std::move(id),
            std::move(rank_sequence),
            internal_id
        );

        ++internal_id;
    }

    if (records.empty()) {
        throw std::runtime_error("The reference file is empty, which is not allowed.\n");
    }

    return references { .records = std::move(records), .total_sequence_length = total_length };
}

queries read_queries(std::filesystem::path const& queries_path) {
    spdlog::info("reading queries from {}", queries_path);

    std::vector<query_record> records{};

    size_t total_length = 0;

    for (auto const record_view : ivio::fastq::reader{{ .input = queries_path }}) {
        std::string const id = internal::extract_record_id(record_view.id);

        if (record_view.seq.empty()) {
            spdlog::warn(
                "The record {} in the query file has an empty sequence and will be skipped.\n",
                id
            );

            continue;
        }

        std::vector<uint8_t> const rank_sequence = internal::chars_to_rank_sequence(record_view.seq);
        std::string const quality(record_view.qual);

        assert(record_view.qual.size() == record_view.seq.size());

        total_length += rank_sequence.size();

        records.emplace_back(
            std::move(id),
            std::move(rank_sequence),
            std::move(quality)
        );
    }

    return queries{ .records = std::move(records), .total_sequence_length = total_length };
}

fmindex load_index(std::filesystem::path const& index_path) {
    spdlog::info("loading index from {}", index_path);

    auto ifs     = std::ifstream(index_path, std::ios::binary);
    auto archive = cereal::BinaryInputArchive{ifs};
    auto index = fmindex{};
    archive(index);

    return index;
}

namespace internal {

std::string extract_record_id(std::string_view const& reference_name) {
    return std::string(reference_name.begin(), std::ranges::find(reference_name, ' '));
}

std::vector<uint8_t> chars_to_rank_sequence(std::string_view const sequence) {
    auto rank_sequence = ivs::convert_char_to_rank<ivs::d_dna5>(sequence);

    uint8_t const replacement_rank = ivs::d_dna5::char_to_rank('N');
    std::ranges::replace_if(
        rank_sequence,
        [] (uint8_t const rank) { return !ivs::verify_rank(rank); },
        replacement_rank
    );

    return rank_sequence;
}

} // namespace internal

} // namespace input
