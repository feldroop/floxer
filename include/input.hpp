#pragma once

#include <floxer_cli.hpp>
#include <fmindex.hpp>

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include <ivio/ivio.h>

namespace input {

struct reference_record {
    std::string const id;
    std::vector<uint8_t> const rank_sequence;
    size_t const internal_id;
};

struct query_record {
    std::string id;
    std::vector<uint8_t> rank_sequence;
    std::string quality;
};

struct references {
    std::vector<reference_record> records;
    size_t total_sequence_length;
};

class queries {
public:
    queries(cli::command_line_input const& cli_input);

    std::optional<query_record> next();

private:
    static constexpr size_t MAX_ALLOWED_QUERY_LENGTH = 100'000;

    ivio::fastq::reader reader;
    cli::command_line_input const& cli_input;
};

references read_references(std::filesystem::path const& reference_sequence_path);

fmindex load_index(std::filesystem::path const& _index_path);

// the number of errors allowed for this a queries alignment (edit distance)
// it was either directly given by the user, or is calculated using the given
// error probability
size_t num_errors_from_user_config(size_t const query_length, cli::command_line_input const& cli_input);

namespace internal {

// it is assumed that the record id is the start of the tag until the first space
std::string extract_record_id(std::string_view const& record_tag);

// convert ASCII DNA chars to a rank sequence of ints from 0 to 5.
// all invalid chars become 5 and sentinel '$' becomes 0.
// this means that this program currently can't accurately handle IUPAC degenerate chars
std::vector<uint8_t> chars_to_rank_sequence(std::string_view const chars);

} // namespace internal

} // namespace input
