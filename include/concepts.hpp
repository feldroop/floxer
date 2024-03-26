#pragma once

#include <ranges>

template <typename T>
concept Sequence = std::ranges::sized_range<T>
    && std::ranges::random_access_range<T>
    && requires(T t) {
        {*t.begin()} -> std::common_with<uint8_t>;
    };


template <typename T>
concept Sequences = std::ranges::sized_range<T>
    && std::ranges::random_access_range<T>
    && requires(T t) {
        {*t.begin()} -> Sequence;
    };
