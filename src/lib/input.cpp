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

#include <ivsigma/ivsigma.h>

#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/std.h>
#include <spdlog/spdlog.h>

namespace input {

using floxer_alphabet_t = ivs::d_dna5;

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

queries::queries(cli::command_line_input const& cli_input_)
    : reader{{ .input = cli_input_.queries_path() }},
    num_queries_read{0},
    cli_input{cli_input_} {}

std::optional<query_record> queries::next() {
    while (true) {
        std::optional const record_view_opt = reader.next();

        if (!record_view_opt.has_value()) {
            return std::nullopt;
        }

        auto const& record_view = *record_view_opt;

        std::string const id = internal::extract_record_id(record_view.id);

        if (record_view.seq.empty()) {
            spdlog::warn(
                "The record {} in the query file has an empty sequence and will be skipped.\n",
                id
            );

            continue;
        }

        size_t const sequence_length = record_view.seq.size();

        if (sequence_length > MAX_ALLOWED_QUERY_LENGTH) {
            spdlog::warn("skipping too large query: {}", id);

            continue;
        }

        // two cases that likely don't occur in practice where the errors are configured in a way such that the
        // alignment algorithm makes no sense and floxer just flags them as unaligned
        size_t const query_num_errors = num_errors_from_user_config(sequence_length, cli_input);
        if (
            sequence_length <= query_num_errors ||
            query_num_errors < cli_input.pex_seed_num_errors()
        ) {
            spdlog::warn(
                "skipping query: {} due to bad configuration regarding the number of errors.\n"
                "\tquery length: {}, errors in query: {}, PEX seed errors: {}",
                id,
                sequence_length,
                query_num_errors,
                cli_input.pex_seed_num_errors()
            );

            continue;
        }

        std::vector<uint8_t> const rank_sequence = internal::chars_to_rank_sequence(record_view.seq);
        std::vector<uint8_t> const reverse_complement_rank_sequence = ivs::reverse_complement_rank<floxer_alphabet_t>(rank_sequence);

        std::string const quality(record_view.qual);

        assert(record_view.qual.size() == sequence_length);

        ++num_queries_read;

        return std::make_optional(query_record {
            .id = std::move(id),
            .rank_sequence = std::move(rank_sequence),
            .reverse_complement_rank_sequence = std::move(reverse_complement_rank_sequence),
            .quality = std::move(quality),
            .internal_id = num_queries_read - 1
        });
    }
}

fmindex load_index(std::filesystem::path const& index_path) {
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
    auto rank_sequence = ivs::convert_char_to_rank<floxer_alphabet_t>(sequence);

    uint8_t const replacement_rank = floxer_alphabet_t::char_to_rank('N');
    std::ranges::replace_if(
        rank_sequence,
        [] (uint8_t const rank) { return !ivs::verify_rank(rank); },
        replacement_rank
    );

    return rank_sequence;
}

} // namespace internal

} // namespace input
