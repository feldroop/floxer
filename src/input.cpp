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
#include <spdlog/spdlog.h>

namespace input {

size_t query_record::num_errors_from_user_config(cli::command_line_input const& cli_input) const {
    return cli_input.query_error_probability().has_value() ?
        static_cast<size_t>(
            std::ceil(rank_sequence.size() * cli_input.query_error_probability().value())
        ) :
        cli_input.query_num_errors().value();
}

references read_references(std::filesystem::path const& reference_sequence_path) {
    std::vector<reference_record> records{};
    
    size_t id = 0;
    size_t num_unnamed = 0;
    size_t total_length = 0;

    std::unordered_map<std::string, size_t> reference_names{};
    bool duplicate_name_warning_given = false;

    for (auto const record_view : ivio::fasta::reader{{ .input = reference_sequence_path }}) {
        std::string raw_tag(record_view.id);
        if (raw_tag.empty()) {
            raw_tag = fmt::format("_reference_without_name_{}_", num_unnamed);
            ++num_unnamed;
        }

        if (record_view.seq.empty()) {
            spdlog::warn(
                "The record {} in the reference file has an empty sequence and will be skipped.\n",
                raw_tag
            );

            continue;
        }

        std::string sam_format_sanitized_name = internal::sanitize_reference_name_for_sam(raw_tag);

        if (reference_names.contains(sam_format_sanitized_name)) {
            if (!duplicate_name_warning_given) {
                spdlog::warn(
                    "Found duplicate names in the reference file. "
                    "These records will be treated separately and given unique names in the output.\n"
                );
                duplicate_name_warning_given = true;
            }
            size_t& num_encountered = reference_names[sam_format_sanitized_name];
            ++num_encountered;
            sam_format_sanitized_name += fmt::format("_{}", num_encountered);
        } else {
            reference_names.emplace(sam_format_sanitized_name, 1);
        }

        std::string const no_degenerate_char_sequence = internal::replace_degenerate_chars(record_view.seq);
        std::vector<uint8_t> const sequence = ivs::convert_char_to_rank<ivs::d_dna5>(no_degenerate_char_sequence);

        auto const result = ivs::verify_rank(sequence);
        if (result.has_value()) {
            size_t const position = result.value();

            throw std::runtime_error(
                fmt::format(
                    "The reference sequence {} "
                    "contians the invalid character {} "
                    "at position {}.\n",
                    record_view.id,
                    record_view.seq[position],
                    position
                )
            );
        }

        spdlog::debug("read reference: \"{}\", length {:L}", raw_tag, sequence.size());

        total_length += sequence.size();
        records.emplace_back(
            id,
            std::move(raw_tag),
            std::move(sam_format_sanitized_name),
            std::move(sequence)
        );

        ++id;
    }

    return references { .records = std::move(records), .total_sequence_length = total_length };
}

queries read_queries(std::filesystem::path const& queries_path) {
    std::vector<query_record> records{};

    size_t id = 0;
    size_t num_unnamed = 0;
    size_t total_length = 0;

    std::unordered_map<std::string, size_t> query_names{};
    bool duplicate_name_warning_given = false;

    for (auto const record_view : ivio::fastq::reader{{ .input = queries_path }}) {
        std::string raw_tag(record_view.id);
        if (raw_tag.empty()) {
            raw_tag = fmt::format("_query_without_name_{}_", num_unnamed);
            ++num_unnamed;
        }

        if (record_view.seq.empty()) {
            spdlog::warn(
                "The record {} in the reference file has an empty sequence and will be skipped.\n",
                raw_tag
            );

            continue;
        }

        std::string sam_format_sanitized_name = internal::sanitize_query_name_for_sam(raw_tag);

        if (query_names.contains(sam_format_sanitized_name)) {
            if (!duplicate_name_warning_given) {
                spdlog::warn(
                    "Found duplicate names in the query file. "
                    "These records will be treated separately and given unique names in the output.\n"
                );
                duplicate_name_warning_given = true;
            }
            size_t& num_encountered = query_names[sam_format_sanitized_name];
            ++num_encountered;
            sam_format_sanitized_name += fmt::format("_{}", num_encountered);
        } else {
            query_names.emplace(sam_format_sanitized_name, 1);
        }

        std::string quality(record_view.qual);
        if (quality.size() != record_view.seq.size()) {
            spdlog::warn(
                "The quality of record {} does not have "
                "the correct length and will be ignored.\n",
                raw_tag
            );
            quality.clear();
        }

        std::string const no_degenerate_char_sequence = internal::replace_degenerate_chars(record_view.seq);
        std::vector<uint8_t> const rank_sequence = ivs::convert_char_to_rank<ivs::d_dna5>(no_degenerate_char_sequence);

        auto const result = ivs::verify_rank(rank_sequence);
        if (result.has_value()) {
            size_t const position = result.value();

            spdlog::warn(
                "Skipped the query {} "
                "due to the invalid character {} "
                "at position {}.\n",
                raw_tag,
                record_view.seq[position],
                position
            );

            continue;
        }

        total_length += rank_sequence.size();
        records.emplace_back(
            id,
            std::move(raw_tag),
            std::move(sam_format_sanitized_name),
            std::move(rank_sequence),
            std::move(quality)
        );

        ++id;
    }

    return queries{ .records = std::move(records), .total_sequence_length = total_length };
}

fmindex load_index(std::filesystem::path const& _index_path) {
    auto ifs     = std::ifstream(_index_path, std::ios::binary);
    auto archive = cereal::BinaryInputArchive{ifs};
    auto index = fmindex{};
    archive(index);
    return index;
}

namespace internal {

std::string sanitize_reference_name_for_sam(std::string const& reference_name) {
    static constexpr char replacement_char = '_';

    std::string sanitized_name = reference_name;
    size_t const valid_start_index = reference_name.find_first_not_of("*=");

    for (size_t i = 0; i < valid_start_index; ++i) {
        sanitized_name[i] = replacement_char;
    }

    std::ranges::replace_if(sanitized_name, [] (char const& old_char) {
        static const std::string invalid_string = "\\,\"\'`()[]{}<>";
        static const std::unordered_set<char> invalid_chars(
            invalid_string.begin(),
            invalid_string.end()
        );
        return old_char < '!' || old_char > '~' || invalid_chars.contains(old_char);
    }, replacement_char);

    return sanitized_name;
}

std::string sanitize_query_name_for_sam(std::string const& query_name) {
    static constexpr char replacement_char = '_';
    static constexpr size_t sam_max_allowed_length = 254;
    
    std::string sanitized_name = query_name.substr(0, sam_max_allowed_length);

    std::ranges::replace_if(sanitized_name, [] (char const& old_char) {
        static constexpr char invalid_char = '@';
        return old_char < '!' || old_char > '~' || old_char == invalid_char;
    }, replacement_char);

    return sanitized_name;
}

constexpr std::array<char, 256> degenerate_to_simple_char_conversion_table() {
    std::array<char, 256> conversion;
    std::iota(conversion.begin(), conversion.end(), 0);

    // arbitrary definition to get rid of degenerate chars for now:
    // convert all degenerate chars except for N to the smallest
    // character (in the order A,C,G,T) that they could be
    // example: d/D could be A,G or T and becomes A
    conversion['r'] = 'a', conversion['R'] = 'A';
    conversion['y'] = 'c', conversion['Y'] = 'C';
    conversion['k'] = 'g', conversion['K'] = 'G';
    conversion['m'] = 'a', conversion['M'] = 'A';
    conversion['s'] = 'c', conversion['S'] = 'C';
    conversion['w'] = 'a', conversion['W'] = 'A';
    conversion['b'] = 'c', conversion['B'] = 'C';
    conversion['d'] = 'a', conversion['D'] = 'A';
    conversion['h'] = 'a', conversion['H'] = 'A';
    conversion['v'] = 'a', conversion['V'] = 'A';

    return conversion;
}

std::string replace_degenerate_chars(std::string_view const& sequence) {
    static constexpr std::array<char, 256> conversion = degenerate_to_simple_char_conversion_table();

    auto const replaced_view = sequence | std::views::transform([&] (char const& c) {
        return conversion[c];
    });

    return std::string(replaced_view.begin(), replaced_view.end());
}

} // namespace internal

} // namespace input
