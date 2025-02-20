#include <about_floxer.hpp>
#include <input.hpp>
#include <search.hpp>

#include <random>
#include <span>
#include <vector>

#include <fmindex-collection/search/SearchNg21.h>
#include <sharg/all.hpp>
#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/ranges.h>

// wrapped in an extra vector for fmindex lib API compatibility
std::vector<std::vector<uint8_t>> create_random_sequence(size_t const length, std::mt19937& random_generator) {
    std::uniform_int_distribution<uint8_t> base_rank_distribution(0, 3);
    std::vector<uint8_t> seq{};

    for (size_t i = 0; i < length; ++i) {
        seq.push_back(base_rank_distribution(random_generator));
    }

    return std::vector<std::vector<uint8_t>>{ std::move(seq) };
}

int main(int argc, char** argv)  {
    sharg::parser parser{ "count_random_occurrences", argc, argv, sharg::update_notifications::off };

    parser.info.author = about_floxer::author;
    parser.info.description = {
        "Search random strings in an FM-Index and output statistics on how often they were found on average."
    };
    parser.info.email = about_floxer::email;
    parser.info.url = about_floxer::url;
    parser.info.short_description = "Search random strings in an FM-Index";
    parser.info.synopsis = {
        "./count_random_occurrences --index index.flxi",
    };
    parser.info.version = "1.0.0";
    parser.info.date = "13.05.2024";

    std::filesystem::path index_path{};
    size_t min_length = 10;
    size_t max_length = 60;
    size_t num_searches_per_length = 1'000'000;
    size_t min_errors = 0;
    size_t max_errors = 3;

    parser.add_option(index_path, sharg::config{
        .short_id = 'i',
        .long_id = "index",
        .description = "The FM-Index file in which to search (created by floxer).",
        .required = true,
        .validator = sharg::input_file_validator{}
    });

    parser.add_option(min_length, sharg::config{
        .short_id = 'm',
        .long_id = "min-length",
        .description = "The shortest length of a random pattern that will be searched."
    });

    parser.add_option(max_length, sharg::config{
        .short_id = 'n',
        .long_id = "max-length",
        .description = "The biggest length of a random pattern that will be searched."
    });

    parser.add_option(num_searches_per_length, sharg::config{
        .short_id = 's',
        .long_id = "searches",
        .description = "The number of searches per length/error combination."
    });

    parser.add_option(min_errors, sharg::config{
        .short_id = 'e',
        .long_id = "min-errors",
        .description = "The minimum nuber of errors with which to search."
    });

    parser.add_option(max_errors, sharg::config{
        .short_id = 'x',
        .long_id = "max-errors",
        .description = "The maximum number with which to search."
    });

    parser.parse();

    auto index = input::load_index(index_path);

    std::mt19937 random_generator(837103474);
    search::internal::search_scheme_cache scheme_cache;

    fmt::print("runs = [\n");

    for (size_t num_errors = min_errors; num_errors <= max_errors; ++num_errors) {
        std::vector<double> count_averages{};

        for (size_t length = min_length; length <= max_length; ++length) {
            double average_count = 0;

            for (size_t i = 0; i < num_searches_per_length; ++i) {
                auto const seq = create_random_sequence(length, random_generator);

                auto const& search_scheme = scheme_cache.get(
                    length,
                    num_errors
                );

                size_t count = 0;

                fmindex_collection::search_ng21::search(
                    index,
                    std::span(seq),
                    search_scheme,
                    [&count] (
                        [[maybe_unused]] size_t const _seed_index_in_wrapper_range,
                        auto cursor,
                        [[maybe_unused]] size_t const _errors
                    ) {
                        count += cursor.count();
                    }
                );

                average_count += count / static_cast<double>(num_searches_per_length);
            }

            count_averages.push_back(average_count);
        }

        fmt::print("    {{ num_errors = {}, count_averages = {} }},\n", num_errors,  count_averages);
    }

    fmt::print("]\n");

    return 0;
}
