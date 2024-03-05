#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace verification {

enum class alignment_variant {
    match, deletion, insertion
};

struct query_alignment {
    // half open range [start_in_reference, end_in_reference)
    size_t start_in_reference;
    size_t end_in_reference;
    size_t num_errors;

    // the terminolgy ov the alignment variants refers is from the perspective of the query
    // e.g. deletion -> deletion in the query
    std::vector<alignment_variant> alignment;

    size_t length_in_reference();
};

std::vector<alignment_variant> alignment_from_string(std::string const& s);

// if at least one adequate alignment is found, the shortest one is returned
// in case of ties, the earliest in the reference is returned
std::optional<query_alignment> query_occurs(
    std::span<const uint8_t> reference,
    std::span<const uint8_t> query,
    size_t const num_allowed_errors
);

} // namespace verification
