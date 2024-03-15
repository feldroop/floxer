#include <input.hpp>

#include <algorithm>
#include <cmath>
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
    return std::isnan(opt.query_error_probability) ?
        opt.query_num_errors : 
        static_cast<size_t>(
            std::ceil(sequence_length * opt.query_error_probability)
        );
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
    
    for (auto const record_view : ivio::fasta::reader{{ .input = reference_sequence_path }}) {
        std::string const raw_tag(record_view.id);
        std::string const sam_format_sanitized_name = sanitize_reference_name_for_sam(raw_tag);
        std::vector<uint8_t> const sequence = ivs::convert_char_to_rank<ivs::d_dna4>(record_view.seq);

        auto const result = ivs::verify_rank(sequence);
        if (result.has_value()) {
            size_t const position = result.value();

            fmt::print(
                stderr, 
                "[INPUT ERROR]\nThe reference sequence {} "
                "contians the invalid character {} "
                "at position {}.\n",
                record_view.id,
                sequence[position],
                position
            );

            exit(-1);
        }

        size_t const sequence_length = sequence.size();
        records.emplace_back(
            std::move(raw_tag),
            std::move(sam_format_sanitized_name),
            std::move(sequence),
            sequence_length
        );
    }

    return records;
}

std::vector<query_record> read_queries(std::filesystem::path const& queries_path) {
    std::vector<query_record> records{};

    for (auto const record_view : ivio::fastq::reader{{ .input = queries_path }}) {
        std::string const raw_tag(record_view.id);
        std::string const sam_format_sanitized_name = sanitize_query_name_for_sam(raw_tag);
        std::string const quality(record_view.qual);
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
                record_view.id,
                char_sequence[position],
                position
            );

            continue;
        }
        size_t const sequence_length = char_sequence.size();
        records.emplace_back(
            std::move(raw_tag),
            std::move(sam_format_sanitized_name),
            std::move(rank_sequence),
            std::move(char_sequence),
            std::move(quality),
            sequence_length
        );
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
