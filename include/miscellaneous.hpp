#pragma once

#include <cstddef>
#include <tuple>
#include <unordered_map>
#include <utility>

// ---------- hash for tuples with hashable members ----------
// inspired by: https://codereview.stackexchange.com/questions/136770/hashing-a-tuple-in-c17

// helper function that does the actual work
template<typename Tuple, size_t... ids>
size_t hash_tuple(Tuple const& tup, std::index_sequence<ids...>) {
    size_t state = 0;
    auto const h = [&]<typename T> (T const& elem) {
        return std::hash<T>()(elem);
    };
    
    for (size_t const hash_value : { h(std::get<ids>(tup))... }) {
        state ^= hash_value * 0x9e3779b9 + (state << 6) + (state >> 1);    
    }

    return state;
}

// hash specialisation
template<typename... Ts>
struct std::hash<std::tuple<Ts...>> {
    size_t operator()(std::tuple<Ts...> const& tup) const noexcept {
        return hash_tuple(tup, std::index_sequence_for<Ts...>());
    }
};

inline size_t ceil_div(size_t const a, size_t const b) {
    return (a % b) ? a / b + 1 : a / b;
}