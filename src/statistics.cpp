#include <statistics.hpp>

namespace statistics {

search_and_alignment_statistics::histogram::histogram(
    histogram_config const& config_,
    std::string const& name_
) : config{config_}, name{name_} {
    data.resize(config.thresholds.size() + 1, 0);
}

void search_and_alignment_statistics::histogram::add_value(size_t const value) {
    ++num_values;

    for (size_t i = 0; i < config.thresholds.size(); ++i) {
        if (value <= config.thresholds[i]) {
            ++data[i];
            return;
        }
    }

    ++data.back();
}

void search_and_alignment_statistics::histogram::merge_with(histogram const& other) {
    assert(config.thresholds == other.donfig.thresholds);
    num_values += other.num_values;

    for (size_t i = 0; i < data.size(); ++i) {
        data[i] += other.data[i];
    }
}

std::string search_and_alignment_statistics::histogram::format_to_string() const {
    return fmt::format(
        "histogram for {} (total: {})\n"
        "threshold:\t{}\tinf\n"
        "occurrences:\t{}",
        name,
        num_values,
        fmt::join(config.thresholds, "\t"),
        fmt::join(data, "\t")
    );
}

void search_and_alignment_statistics::insert_value_to(
    std::string const& target_name,
    size_t const value
) {
    auto iter = std::ranges::find_if(histograms, [&] (histogram const& histo) {
        return histo.name == target_name;
    });

    if (iter == histograms.end()) {
        throw std::runtime_error("Internal bug in stats generation");
    }

    iter->add_value(value);
}

void search_and_alignment_statistics::add_query_length(size_t const value) {
    insert_value_to(query_lengths_name, value);
}

void search_and_alignment_statistics::add_seed_length(size_t const value) {
    insert_value_to(seed_lengths_name, value);
}

void search_and_alignment_statistics::add_num_anchors_per_seed(size_t const value) {
    insert_value_to(anchors_per_seed_name, value);
}

void search_and_alignment_statistics::add_num_anchors_per_query(size_t const value) {
    insert_value_to(anchors_per_query_name, value);
}

void search_and_alignment_statistics::add_num_alignments(size_t const value) {
    insert_value_to(alignments_per_query_name, value);
}

void search_and_alignment_statistics::add_alignment_edit_distance(size_t const value) {
    insert_value_to(alignments_edit_distance_name, value);
}

std::vector<std::string> search_and_alignment_statistics::format_histograms() const {
    std::vector<std::string> formatted_histograms{};
    
    for (auto const& histo : histograms) {
        formatted_histograms.push_back(histo.format_to_string());
    }

    return formatted_histograms;
}

search_and_alignment_statistics& combine_stats(
    search_and_alignment_statistics & inout,
    search_and_alignment_statistics const& other
) {
    for (size_t i = 0; i < inout.histograms.size(); ++i) {
        inout.histograms.at(i).merge_with(other.histograms.at(i));
    }

    return inout;
}

} // namespace statistics