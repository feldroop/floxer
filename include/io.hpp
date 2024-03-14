#pragma once

#include <about_floxer.hpp>
#include <fmindex.hpp>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace io {

struct reference_record {
    std::string const raw_tag;
    std::string const sam_format_sanitized_name;
    std::vector<uint8_t> const sequence;
    size_t const sequence_length;
};

struct query_record {
    std::string const raw_tag;
    std::string const sam_format_sanitized_name;
    std::vector<uint8_t> const sequence;
    std::string const quality;
    size_t const sequence_length;
};

std::vector<reference_record> read_references(std::filesystem::path const& reference_sequence_path);

std::vector<query_record> read_queries(std::filesystem::path const& queries_path);

void save_index(fmindex const& _index, std::filesystem::path const& _index_path);

fmindex load_index(std::filesystem::path const& _index_path);

// ---------- sam output ----------
struct file_level_metadata_t {
    std::string const version = "1.6"; // VN
    // no alignement sorting order for now, because floxer only does grouping
    std::string const alignment_grouping = "query"; // GO

    template<class Stream>
    void append_to_stream(Stream& stream) const {
        stream << "@HG\tVN:" << version
            << "\tGO:" << alignment_grouping << '\n';
    }
};

struct reference_sequence_metadata {
    std::string const name; // SN
    uint32_t const length; // LN

    template<class Stream>
    void append_to_stream(Stream& stream) const {
        stream << "@SQ\tSN:" << name
            << "\tLN:" << length << '\n';
    }
};

struct program_t {
    size_t const id = 0; // ID
    std::string const name = "floxer"; // PN
    // no command line call, because it might leak unexpected things
    std::string const description = about_floxer::long_description + " " + about_floxer::url; // DS
    std::string const version = about_floxer::version; // VN

    template<class Stream>
    void append_to_stream(Stream& stream) const {
        stream << "@PG\tID:" << id
            << "\tPN:" << name
            << "\tDS:" << description
            << "\tVN:" << version << '\n';
    }
};

class sam_header {
    file_level_metadata_t const file_level_metadata{}; // @HD
    std::vector<reference_sequence_metadata> reference_sequences{}; // @SQ
    // no @RG read groups for now
    program_t const program; // @PG

public:
    sam_header(std::vector<reference_record> const& reference_records);

    std::string to_string() const;
};

} // namespace io
