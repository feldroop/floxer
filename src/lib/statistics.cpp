#include <statistics.hpp>

#include <algorithm>

#include <spdlog/fmt/ranges.h>

namespace statistics {

std::string search_and_alignment_statistics::count::format_to_string_for_stdout() const {
    return fmt::format("number of {}: {}", name, value);
}

std::string search_and_alignment_statistics::count::format_as_toml() const {
    std::string underscored_name = name;
    std::ranges::replace(underscored_name, ' ', '_');

    return fmt::format("{} = {}\n", underscored_name, value);
}

search_and_alignment_statistics::histogram::histogram(
    histogram_config const& config_,
    std::string const& name_
) : config{config_}, name{name_} {
    data.resize(config.thresholds.size() + 1, 0);
}

void search_and_alignment_statistics::histogram::add_value(size_t const value) {
    ++num_values;
    min = std::min(min, value);
    sum += value;
    max = std::max(max, value);

    for (size_t i = 0; i < config.thresholds.size(); ++i) {
        if (value <= config.thresholds[i]) {
            ++data[i];
            return;
        }
    }

    ++data.back();
}

void search_and_alignment_statistics::histogram::merge_with(histogram const& other) {
    assert(config.thresholds == other.config.thresholds);
    num_values += other.num_values;
    min = std::min(min, other.min);
    sum += other.sum;
    max = std::max(max, other.max);

    for (size_t i = 0; i < data.size(); ++i) {
        data[i] += other.data[i];
    }
}

std::string search_and_alignment_statistics::histogram::format_to_string_for_stdout() const {
    std::string basic_stats = num_values > 0 ?
        fmt::format(
            "\nmin = {}, mean = {:.2f}, max = {}", min, sum / num_values, max
        ) : "";

    return fmt::format(
        "histogram for {} (total: {})\n"
        "threshold:\t{}\tinf\n"
        "occurrences:\t{}"
        "{}",
        name,
        num_values,
        fmt::join(config.thresholds, "\t"),
        fmt::join(data, "\t"),
        basic_stats
    );
}

std::string search_and_alignment_statistics::histogram::format_as_toml() const {
    std::string underscored_name = name;
    std::ranges::replace(underscored_name, ' ', '_');

    std::string formatted = fmt::format(
        "[{}]\n"
        "num_values = {}\n"
        "thresholds = {}\n"
        "occurrences = {}\n",
        underscored_name,
        num_values,
        config.thresholds,
        data
    );

    if (num_values > 0) {
        formatted += fmt::format(
            "min_value = {}\n"
            "mean = {:.2f}\n"
            "max_value = {}\n",
            min,
            sum / num_values,
            max
        );
    }

    return formatted;
}

search_and_alignment_statistics::count&
search_and_alignment_statistics::count_by_name(std::string const& name) {
    auto iter = std::ranges::find_if(counts, [&] (count const& c) {
        return c.name == name;
    });

    if (iter == counts.end()) {
        throw std::runtime_error("Internal bug in stats generation");
    }

    return *iter;
}

search_and_alignment_statistics::count const&
search_and_alignment_statistics::count_by_name(std::string const& name) const {
    auto iter = std::ranges::find_if(counts, [&] (count const& c) {
        return c.name == name;
    });

    if (iter == counts.end()) {
        throw std::runtime_error("Internal bug in stats generation");
    }

    return *iter;
}

search_and_alignment_statistics::histogram&
search_and_alignment_statistics::histogram_by_name(std::string const& name) {
    auto iter = std::ranges::find_if(histograms, [&] (histogram const& histo) {
        return histo.name == name;
    });

    if (iter == histograms.end()) {
        throw std::runtime_error("Internal bug in stats generation");
    }

    return *iter;
}

search_and_alignment_statistics::histogram const&
search_and_alignment_statistics::histogram_by_name(std::string const& name) const {
    auto iter = std::ranges::find_if(histograms, [&] (histogram const& histo) {
        return histo.name == name;
    });

    if (iter == histograms.end()) {
        throw std::runtime_error("Internal bug in stats generation");
    }

    return *iter;
}

void search_and_alignment_statistics::increment_count(std::string const& target_name) {
    ++count_by_name(target_name).value;
}

void search_and_alignment_statistics::insert_value_to(
    std::string const& target_name,
    size_t const value
) {
    histogram_by_name(target_name).add_value(value);
}

void search_and_alignment_statistics::increment_num_completely_excluded_queries() {
    increment_count(num_completely_excluded_queries_name);
}

void search_and_alignment_statistics::add_query_length(size_t const value) {
    insert_value_to(query_lengths_name, value);
}

void search_and_alignment_statistics::add_seed_length(size_t const value) {
    insert_value_to(seed_lengths_name, value);
}

void search_and_alignment_statistics::add_num_errors_per_seed(size_t const value) {
    insert_value_to(errors_per_seed_name, value);
}

void search_and_alignment_statistics::add_num_seeds_per_query(size_t const value) {
    insert_value_to(seeds_per_query_name, value);
}

void search_and_alignment_statistics::add_statistics_for_seeds(std::vector<search::seed> const& seeds) {
    add_num_seeds_per_query(seeds.size());

    for (auto const& seed : seeds) {
        add_num_errors_per_seed(seed.num_errors);
        add_seed_length(seed.sequence.size());
    }
}

void search_and_alignment_statistics::add_num_anchors_per_seed(size_t const value) {
    insert_value_to(anchors_per_seed_name, value);
}

void search_and_alignment_statistics::add_num_kept_anchors_per_partly_excluded_seed(size_t const value) {
    insert_value_to(kept_anchors_per_partly_excluded_seed_name, value);
}

void search_and_alignment_statistics::add_num_raw_anchors_per_excluded_seed(size_t const value) {
    insert_value_to(raw_anchors_per_excluded_seed_name, value);
}

void search_and_alignment_statistics::add_num_anchors_per_query(size_t const value) {
    insert_value_to(anchors_per_query_name, value);
}

void search_and_alignment_statistics::add_num_excluded_raw_anchors_per_query(size_t const value) {
    insert_value_to(excluded_raw_anchors_per_query_name, value);
}

void search_and_alignment_statistics::add_reference_span_size_aligned_inner_node(size_t const value) {
    insert_value_to(reference_span_sizes_aligned_inner_nodes_name, value);
}

void search_and_alignment_statistics::add_reference_span_size_avoided_inner_node(size_t const value) {
    insert_value_to(reference_span_sizes_avoided_inner_nodes_name, value);
}

void search_and_alignment_statistics::add_reference_span_size_aligned_root(size_t const value) {
    insert_value_to(reference_span_sizes_aligned_root_name, value);
}

void search_and_alignment_statistics::add_reference_span_size_avoided_root(size_t const value) {
    insert_value_to(reference_span_sizes_avoided_root_name, value);
}

void search_and_alignment_statistics::add_num_alignments(size_t const value) {
    insert_value_to(alignments_per_query_name, value);
}

void search_and_alignment_statistics::add_alignment_edit_distance(size_t const value) {
    insert_value_to(alignments_edit_distance_name, value);
}

void search_and_alignment_statistics::add_statistics_for_search_result(
    search::search_result const& search_result
) {
    size_t num_anchors_of_whole_query = 0;
    size_t num_excluded_anchors_of_whole_query = 0;
    bool all_excluded = true;

    for (auto const& anchors_of_seed : search_result.anchors_by_seed) {
        switch (anchors_of_seed.status) {
            case search::seed_status::fully_excluded :
                add_num_raw_anchors_per_excluded_seed(anchors_of_seed.num_excluded_raw_anchors);
                num_excluded_anchors_of_whole_query += anchors_of_seed.num_excluded_raw_anchors;
                break;
            case search::seed_status::partly_excluded :
                add_num_kept_anchors_per_partly_excluded_seed(anchors_of_seed.num_kept_useful_anchors);
                num_excluded_anchors_of_whole_query += anchors_of_seed.num_excluded_raw_anchors;
                num_anchors_of_whole_query += anchors_of_seed.num_kept_useful_anchors;
                all_excluded = false;
                break;
            case search::seed_status::not_excluded :
                add_num_anchors_per_seed(anchors_of_seed.num_kept_useful_anchors);
                num_anchors_of_whole_query += anchors_of_seed.num_kept_useful_anchors;
                all_excluded = false;
                break;
            default:
                throw std::runtime_error("(should be unreachable) internal bug in statistics gathering");
        }
    }

    add_num_anchors_per_query(num_anchors_of_whole_query);
    add_num_excluded_raw_anchors_per_query(num_excluded_anchors_of_whole_query);
    if (all_excluded) {
        increment_num_completely_excluded_queries();
    }
}

size_t search_and_alignment_statistics::num_queries() const {
    return histogram_by_name(query_lengths_name).num_values;
}

std::vector<std::string> search_and_alignment_statistics::format_statistics_for_stdout() const {
    std::vector<std::string> formatted_statistics{};

    for (auto const& c : counts) {
        formatted_statistics.emplace_back(c.format_to_string_for_stdout());
    }

    for (auto const& histo : histograms) {
        formatted_statistics.emplace_back(histo.format_to_string_for_stdout());
    }

    return formatted_statistics;
}

std::string search_and_alignment_statistics::format_statistics_as_toml() const {
    std::string formatted{};

    for (auto const& c : counts) {
        formatted += c.format_as_toml();
    }

    for (auto const& histo : histograms) {
        formatted += histo.format_as_toml();
    }

    return formatted;
}

search_and_alignment_statistics& combine_stats(
    search_and_alignment_statistics & inout,
    search_and_alignment_statistics const& other
) {
    for (size_t i = 0; i < inout.counts.size(); ++i) {
        inout.counts.at(i).value += other.counts.at(i).value;
    }

    for (size_t i = 0; i < inout.histograms.size(); ++i) {
        inout.histograms.at(i).merge_with(other.histograms.at(i));
    }

    return inout;
}

} // namespace statistics
