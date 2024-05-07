#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>

namespace math {

inline int32_t saturate_value_to_int32_max(size_t const value) {
    if (value > std::numeric_limits<int32_t>::max()) {
        return std::numeric_limits<int32_t>::max();
    } else {
        return static_cast<int32_t>(value);
    }
}

inline size_t ceil_div(size_t const a, size_t const b) {
    return (a % b) ? a / b + 1 : a / b;
}

} // namespace math
