#pragma once

#include <cmath>
#include <filesystem>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>

#include <spdlog/fmt/fmt.h>

namespace cli {

class command_line_input {
    template<typename T>
    struct cli_option {
        const char short_id;
        const std::string long_id;
        T value;

        std::string command_line_call() const {
            if constexpr (std::is_same<T, std::filesystem::path>::value) {
                return fmt::format(
                    " --{} {}{}",
                    long_id,
                    value.has_parent_path() ? ".../" : "",
                    value.filename().c_str()
                );
            } else if constexpr (std::is_same<T, bool>::value) {
                return fmt::format(" --{}", long_id);
            } else {
                return fmt::format(" --{} {}", long_id, value);
            }
        }
    };

    cli_option<std::filesystem::path> reference_path_{ 'r', "reference", "" };
    cli_option<std::filesystem::path> queries_path_{ 'q', "queries", "" };
    cli_option<std::filesystem::path> output_path_{ 'o', "output", "" };
    cli_option<std::filesystem::path> index_path_{ 'i', "index", "" };
    cli_option<std::filesystem::path> logfile_path_{ 'l', "logfile", "" };

    cli_option<size_t> query_num_errors_{ 'e', "query-errors", std::numeric_limits<size_t>::max() };
    cli_option<double> query_error_probability_{ 'p', "error-probability", NAN };
    cli_option<size_t> pex_seed_num_errors_{ 's', "seed-errors", 2 };
    cli_option<size_t> num_threads_{ 't', "threads", 1 };
    
    cli_option<bool> print_stats_{ 'a', "print-stats", false };

public:
    void parse_and_validate(int argc, char ** argv); 
    
    std::filesystem::path const& reference_path() const;
    std::filesystem::path const& queries_path() const;
    std::filesystem::path const& output_path() const;
    std::optional<std::filesystem::path> index_path() const;
    std::optional<std::filesystem::path> logfile_path() const;

    std::optional<size_t> query_num_errors() const;
    std::optional<double> query_error_probability() const;
    size_t pex_seed_num_errors() const;
    size_t num_threads() const;
    bool print_stats() const;

    // not the exact one, but a sanitized and canonical version
    std::string command_line_call() const;

private:
    void validate() const;
};

} // namespace cli
