#include <input.hpp>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <unordered_set>

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>

#include <fmt/core.h>

namespace input {

size_t query_record::num_errors_from_user_config(cli::options const& opt) const {
    return opt.query_error_probability_was_set() ?
        static_cast<size_t>(
            std::ceil(sequence_length * opt.query_error_probability)
        ) :
        opt.query_num_errors;
}

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

std::vector<reference_record> read_references(std::filesystem::path const& reference_sequence_path) {
    std::vector<reference_record> records{};
    
    size_t id = 0;
    size_t num_unnamed = 0;

    std::unordered_map<std::string, size_t> reference_names{};
    bool duplicate_name_warning_given = false;

    for (auto const record_view : ivio::fasta::reader{{ .input = reference_sequence_path }}) {
        std::string raw_tag(record_view.id);
        if (raw_tag.empty()) {
            raw_tag = fmt::format("_reference_without_name_{}_", num_unnamed);
            ++num_unnamed;
        }

        if (record_view.seq.empty()) {
            fmt::print(
                stderr,
                "[INPUT WARNING]\nThe record {} in the reference file has an empty sequence and will be skipped.\n",
                raw_tag
            );

            continue;
        }

        std::string sam_format_sanitized_name = sanitize_reference_name_for_sam(raw_tag);

        if (reference_names.contains(sam_format_sanitized_name)) {
            if (!duplicate_name_warning_given) {
                fmt::print(
                    stderr,
                    "[INPUT WARNING]\nFound duplicate names in the reference file. "
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

        std::vector<uint8_t> const sequence = ivs::convert_char_to_rank<ivs::d_dna4>(record_view.seq);

        auto const result = ivs::verify_rank(sequence);
        if (result.has_value()) {
            size_t const position = result.value();

            throw std::runtime_error(
                fmt::format(
                    "[INPUT ERROR]\nThe reference sequence {} "
                    "contians the invalid character {} "
                    "at position {}.\n",
                    record_view.id,
                    sequence[position],
                    position
                )
            );
        }

        size_t const sequence_length = sequence.size();
        records.emplace_back(
            id,
            std::move(raw_tag),
            std::move(sam_format_sanitized_name),
            std::move(sequence),
            sequence_length
        );

        ++id;
    }

    return records;
}

std::vector<query_record> read_queries(std::filesystem::path const& queries_path) {
    std::vector<query_record> records{};

    size_t id = 0;
    size_t num_unnamed = 0;

    std::unordered_map<std::string, size_t> query_names{};
    bool duplicate_name_warning_given = false;

    for (auto const record_view : ivio::fastq::reader{{ .input = queries_path }}) {
        std::string raw_tag(record_view.id);
        if (raw_tag.empty()) {
            raw_tag = fmt::format("_query_without_name_{}_", num_unnamed);
            ++num_unnamed;
        }

        if (record_view.seq.empty()) {
            fmt::print(
                stderr,
                "[INPUT WARNING]\nThe record {} in the reference file has an empty sequence and will be skipped.\n",
                raw_tag
            );

            continue;
        }

        std::string sam_format_sanitized_name = sanitize_query_name_for_sam(raw_tag);

        if (query_names.contains(sam_format_sanitized_name)) {
            if (!duplicate_name_warning_given) {
                fmt::print(
                    stderr,
                    "[INPUT WARNING]\nFound duplicate names in the query file. "
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
            fmt::print(
                stderr,
                "[INPUT WARNING]\nThe quality of record {} does not have "
                "the correct length and will be ignored.\n",
                raw_tag
            );
            quality.clear();
        }

        std::string const char_sequence(record_view.seq);
        std::vector<uint8_t> const rank_sequence = ivs::convert_char_to_rank<ivs::d_dna4>(record_view.seq);

        auto const result = ivs::verify_rank(rank_sequence);
        if (result.has_value()) {
            size_t const position = result.value();

            fmt::print(
                stderr,
                "[INPUT WARNING]\nSkipped the query {} "
                "due to the invalid character {} "
                "at position {}.\n",
                raw_tag,
                char_sequence[position],
                position
            );

            continue;
        }
        size_t const sequence_length = char_sequence.size();
        records.emplace_back(
            id,
            std::move(raw_tag),
            std::move(sam_format_sanitized_name),
            std::move(rank_sequence),
            std::move(char_sequence),
            std::move(quality),
            sequence_length
        );

        ++id;
    }

    return records;
}

fmindex load_index(std::filesystem::path const& _index_path) {
    auto ifs     = std::ifstream(_index_path, std::ios::binary);
    auto archive = cereal::BinaryInputArchive{ifs};
    auto index = fmindex{};
    archive(index);
    return index;
}

} // namespace input
