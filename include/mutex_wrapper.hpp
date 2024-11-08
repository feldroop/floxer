#pragma once

#include <mutex>
#include <shared_mutex>
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

template<typename T>
class shared_mutex_guarded {
public:
    struct write_lock {
        std::unique_lock<std::shared_mutex> guard;
        T& ref;
    };

    struct read_lock {
        std::shared_lock<std::shared_mutex> guard;
        T const& ref;
    };

    shared_mutex_guarded() = default;

    template<class... Args>
    shared_mutex_guarded(Args&&... args) : inner(std::forward<Args>(args)...) {}

    write_lock lock_unique() {
        return write_lock {
            .guard = std::unique_lock<std::shared_mutex>(mutex),
            .ref = inner
        };
    }

    read_lock lock_shared() const {
        return read_lock {
            .guard = std::shared_lock<std::shared_mutex>(mutex),
            .ref = inner
        };
    }

private:
    T inner;
    mutable std::shared_mutex mutex;
};
