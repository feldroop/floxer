#pragma once

#include <mutex>
#include <utility>

template<typename T>
class mutex_guarded {
public:
    struct lock {
        std::lock_guard<std::mutex> guard;
        T& ref;
    };

    mutex_guarded() = default;

    template<class... Args>
    mutex_guarded(Args&&... args) : inner(std::forward<Args>(args)...) {}

    lock lock_unique() {
        return lock {
            .guard = std::lock_guard<std::mutex>(mutex),
            .ref = inner
        };
    }

private:
    T inner;
    mutable std::mutex mutex;
};
