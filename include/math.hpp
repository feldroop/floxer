#pragma once

#include <cmath>
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

inline size_t floating_point_error_aware_ceil(double const value) {
    static constexpr double epsilon = 0.000000001;
    // subtract epsilson such that the ceiling doesn't add 1 for something like 5.000000001
    // add epsilon such that the truncation cast to size_t doesn't subtract 1 for something like 4.999999998
    return std::ceil(value - epsilon) + epsilon;
}

} // namespace math
