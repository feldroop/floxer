#include <input.hpp>
#include <math.hpp>

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

size_t num_errors_from_user_config(size_t const query_length, cli::command_line_input const& cli_input) {
    if (cli_input.query_error_probability().has_value()) {
        return math::floating_point_error_aware_ceil(
            query_length * cli_input.query_error_probability().value()
        );
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

queries read_queries(cli::command_line_input const& cli_input) {
    static constexpr size_t MAX_ALLOWED_QUERY_LENGTH = 100'000;

    spdlog::info("reading queries from {}", cli_input.queries_path());

    std::vector<query_record> records{};
    std::vector<query_record> records_with_invalid_config{};

    size_t total_length = 0;

    for (auto const record_view : ivio::fastq::reader{{ .input = cli_input.queries_path() }}) {
        std::string const id = internal::extract_record_id(record_view.id);

        if (record_view.seq.empty()) {
            spdlog::warn(
                "The record {} in the query file has an empty sequence and will be skipped.\n",
                id
            );

            continue;
        }

        if (record_view.seq.size() > MAX_ALLOWED_QUERY_LENGTH) {
            spdlog::warn("skipping too large query: {}", id);

            continue;
        }

        std::vector<uint8_t> const rank_sequence = internal::chars_to_rank_sequence(record_view.seq);
        std::string const quality(record_view.qual);

        assert(record_view.qual.size() == record_view.seq.size());

        auto record = query_record {
            .id = std::move(id),
            .rank_sequence = std::move(rank_sequence),
            .quality = std::move(quality)
        };

        // two cases that likely don't occur in practice where the errors are configured in a way such that the
        // alignment algorithm makes no sense and floxer just flags them as unaligned
        size_t const query_num_errors = num_errors_from_user_config(record.rank_sequence.size(), cli_input);
        if (
            record.rank_sequence.size() <= query_num_errors ||
            query_num_errors < cli_input.pex_seed_num_errors()
        ) {
            spdlog::debug(
                "skipping query: {} due to bad configuration regarding the number of errors.\n"
                "\tquery length: {}, errors in query: {}, PEX seed errors: {}",
                record.id,
                record.rank_sequence.size(),
                query_num_errors,
                cli_input.pex_seed_num_errors()
            );

            records_with_invalid_config.emplace_back(std::move(record));

            continue;
        }

        total_length += record.rank_sequence.size();
        records.emplace_back(std::move(record));
    }

    return queries{
        .records = std::move(records),
        .records_with_invalid_config = std::move(records_with_invalid_config),
        .total_sequence_length = total_length
    };
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

std::string extract_record_id(std::string_view const& record_tag) {
    return std::string(record_tag.begin(), std::ranges::find(record_tag, ' '));
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
