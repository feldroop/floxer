#include <about_floxer.hpp>

#include <algorithm>
#include <concepts>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

#include <sharg/all.hpp>
#include <spdlog/fmt/fmt.h>
#include <spdlog/spdlog.h>

using namespace seqan3::literals;

// this needs to be added, because the tp tag is used defined by minimap
template <>
struct seqan3::sam_tag_type<"tp"_tag> {
    using type = char;
};

size_t get_max_edit_distance(size_t const sequence_length, double const error_rate) {
    double const num_errors_frac = sequence_length * error_rate;

    // handle floating point inaccuracy
    static constexpr double epsilon = 0.000000001;
    if (std::abs(num_errors_frac - std::round(num_errors_frac)) < epsilon) {
        return static_cast<size_t>(std::round(num_errors_frac) + epsilon);
    } else {
        return static_cast<size_t>(std::ceil(num_errors_frac));
    }
}

struct alignment_record_t {
    size_t const num_query_bases_consumed_by_cigar;
    size_t const num_unclipped_query_bases_consumed_by_cigar;
    size_t const num_reference_bases_consumed_by_cigar;

    size_t const num_soft_clipped_bases;
    size_t const num_hard_clipped_bases;

    size_t const query_length_without_clipped_bases;

    size_t const edit_distance;
    double const edit_distance_error_rate;

    bool const is_inversion;

    size_t const longest_indel;

    // possible extensions:
    // - edit distance including clipped
    // - distance function with affine/convex gap costs
    // - reference id, position

    bool is_soft_clipped() const {
        return num_soft_clipped_bases > 0;
    }

    bool is_hard_clipped() const {
        return num_hard_clipped_bases > 0;
    }

    bool is_clipped() const {
        return is_soft_clipped() || is_hard_clipped();
    }

    size_t num_clipped_bases() const {
        return num_soft_clipped_bases + num_hard_clipped_bases;
    }

    bool is_significantly_clipped(double const floxer_allowed_error_rate) const {
        int64_t const errors_left = get_max_edit_distance(
            num_unclipped_query_bases_consumed_by_cigar,
            floxer_allowed_error_rate
        ) - static_cast<int64_t>(edit_distance);
        return !is_high_edit_distance(floxer_allowed_error_rate) &&
            static_cast<int64_t>(num_clipped_bases() * 0.75) > errors_left;
    }

    bool is_high_edit_distance(double const floxer_allowed_error_rate) const {
        size_t const max_edit_distance = get_max_edit_distance(query_length_without_clipped_bases, floxer_allowed_error_rate);
        return edit_distance > max_edit_distance;
    }
};

struct alignment_data_for_query_t {
    bool is_mapped{false};
    bool is_explicitly_unmapped{false};

    std::optional<alignment_record_t> primary_alignment{};

    std::vector<alignment_record_t> supplementary_alignments{};

    // secondary alignments could actually also be representative of secondary chimeric alignment (and thus not really linear)
    std::vector<alignment_record_t> secondary_linear_basic_alignments{};
    std::vector<alignment_record_t> secondary_linear_high_edit_distance_alignments{};
    std::vector<alignment_record_t> secondary_linear_clipped_alignments{};
    std::vector<alignment_record_t> secondary_inverted_alignments{};
    std::vector<alignment_record_t> secondary_supplementary_alignments{};

    void check_floxer_expectations() const {
        if (!secondary_linear_high_edit_distance_alignments.empty()) {
            spdlog::warn("Unexpected high edit distance alignment in floxer.");
        }

        if (!supplementary_alignments.empty()) {
            spdlog::warn("Unexpected primary supplementary alignment in floxer.");
        }

        if (!secondary_inverted_alignments.empty()) {
            spdlog::warn("Unexpected inverted alignment in floxer.");
        }

        if (!secondary_supplementary_alignments.empty()) {
            spdlog::warn("Unexpected secondary supplementary alignment in floxer.");
        }
    }

    void check_general_expectations(
        std::string_view const query_id,
        std::string_view const aligner_name,
        size_t const full_query_length
    ) const {
        if (!(
            (is_mapped && primary_alignment.has_value()) ||
            is_explicitly_unmapped
        )) {
            spdlog::warn("Inconsistent mapping status in {} alignment of query {}.", aligner_name, query_id);
        }

        visit_all_alignment_records([&] (alignment_record_t const& record) {
            if (full_query_length != record.num_query_bases_consumed_by_cigar + record.num_hard_clipped_bases) {
                spdlog::warn(
                    "Inconsistent query lengths in {} alignment of query {}. "
                    "Query length: {}, CIGAR consumed: {}, hard clipped bases: {}",
                    aligner_name,
                    query_id,
                    full_query_length,
                    record.num_query_bases_consumed_by_cigar,
                    record.num_hard_clipped_bases
                );
            }
        });

        if (!secondary_supplementary_alignments.empty() && !is_multiple_mapping()) {
            spdlog::warn(
                "Unexpected {} secondary supplementary alignment without multiple mapping for query {}.",
                aligner_name, query_id
            );
        }
    }

    void visit_all_alignment_records(std::function<void(alignment_record_t const&)> f) const {
        if (primary_alignment.has_value()) {
            f(primary_alignment.value());
        }

        for (auto const& record : secondary_linear_basic_alignments) {
            f(record);
        }

        for (auto const& record : secondary_linear_clipped_alignments) {
            f(record);
        }

        for (auto const& record : secondary_linear_high_edit_distance_alignments) {
            f(record);
        }

        for (auto const& record : supplementary_alignments) {
            f(record);
        }

        for (auto const& record : secondary_inverted_alignments) {
            f(record);
        }
    }

    std::optional<double> primary_error_rate() const {
        return primary_alignment.has_value() ?
            std::make_optional(primary_alignment.value().edit_distance_error_rate) :
            std::nullopt;
    }

    bool is_multiple_mapping() const {
        return is_mapped &&
            (
                !secondary_linear_basic_alignments.empty() ||
                !secondary_linear_clipped_alignments.empty() ||
                !secondary_linear_high_edit_distance_alignments.empty() ||
                !secondary_inverted_alignments.empty()
            );
    }

    bool has_primary_chimeric() const {
        return is_mapped && !supplementary_alignments.empty();
    }

    bool has_primary_inversion() const {
        return is_mapped && primary_alignment.value().is_inversion;
    }

    bool has_primary_linear() const {
        return is_mapped && !has_primary_chimeric() && !has_primary_inversion();
    }

    bool has_primary_linear_basic(double const floxer_allowed_error_rate) const {
        return has_primary_linear()
            && !primary_alignment.value().is_high_edit_distance(floxer_allowed_error_rate)
            && !primary_alignment.value().is_significantly_clipped(floxer_allowed_error_rate);
    }

    bool best_is_chimeric_or_inversion(double const floxer_allowed_error_rate) const {
        return is_mapped
            && !has_basic(floxer_allowed_error_rate)
            && !best_is_high_edit_distance(floxer_allowed_error_rate)
            && !best_is_significantly_clipped(floxer_allowed_error_rate);
    }

    // significantly clipped is defined as "better" than chimeric or inversion
    bool best_is_significantly_clipped(double const floxer_allowed_error_rate) const {
        if (has_basic(floxer_allowed_error_rate) || best_is_high_edit_distance(floxer_allowed_error_rate)) {
            return false;
        }

        return has_primary_linear() || !secondary_linear_clipped_alignments.empty();
    }

    // high edit distance is defined as "better" than significantly clipped
    bool best_is_high_edit_distance(double const floxer_allowed_error_rate) const {
        if (has_basic(floxer_allowed_error_rate)) {
            return false;
        }

        return (
            has_primary_linear()
            && primary_alignment.value().is_high_edit_distance(floxer_allowed_error_rate)
            && !primary_alignment.value().is_significantly_clipped(floxer_allowed_error_rate)
        ) || !secondary_linear_high_edit_distance_alignments.empty();
    }

    // basic = linear, not high edit distance, not signficantly clipped
    bool has_basic(double const floxer_allowed_error_rate) const {
        return has_primary_linear_basic(floxer_allowed_error_rate)
            || !secondary_linear_basic_alignments.empty();
    }

    // gives values for primary or average amongst secondary basics
    double basic_error_rate(double const floxer_allowed_error_rate) const {
        if (!has_basic(floxer_allowed_error_rate)) {
            throw std::runtime_error("alignment data has no basic alignment");
        }

        if (has_primary_linear_basic(floxer_allowed_error_rate)) {
            return primary_alignment.value().edit_distance_error_rate;
        }

        double sum = 0.0;
        for (auto const& alignment : secondary_linear_basic_alignments) {
            sum += alignment.edit_distance_error_rate;
        }

        return sum / secondary_linear_basic_alignments.size();
    }

    double basic_longest_indel(double const floxer_allowed_error_rate) const {
        if (!has_basic(floxer_allowed_error_rate)) {
            throw std::runtime_error("alignment data has no basic alignment");
        }

        if (has_primary_linear_basic(floxer_allowed_error_rate)) {
            return primary_alignment.value().longest_indel;
        }

        double sum = 0.0;
        for (auto const& alignment : secondary_linear_basic_alignments) {
            sum += alignment.longest_indel;
        }

        return sum / secondary_linear_basic_alignments.size();
    }
};

struct query_data_t {
    std::optional<seqan3::dna5_vector> sequence{};

    bool mentioned_by_floxer{false};
    bool mentioned_by_minimap{false};

    alignment_data_for_query_t floxer_alignments{};
    alignment_data_for_query_t minimap_alignments{};

    void check_expectations(std::string_view const query_id) const {
        if (!mentioned_by_floxer) {
            spdlog::warn("Query {} not mentioned by floxer", query_id);
        }

        if (!mentioned_by_minimap) {
            spdlog::warn("Query {} not mentioned by minimap", query_id);
        }

        if (!sequence.has_value()) {
            spdlog::warn("Query {} did not contain full sequence in any aligner file.", query_id);
        }

        floxer_alignments.check_floxer_expectations();
        floxer_alignments.check_general_expectations(query_id, "floxer", sequence.value().size());

        minimap_alignments.check_general_expectations(query_id, "minimap", sequence.value().size());
    }

    bool is_unmapped_floxer() const {
        return floxer_alignments.is_explicitly_unmapped;
    }

    bool is_unmapped_minimap() const {
        return minimap_alignments.is_explicitly_unmapped;
    }

    bool is_unmapped_both() const {
        return is_unmapped_floxer() && is_unmapped_minimap();
    }

    bool is_mapped_both() const {
        return !is_unmapped_floxer() && !is_unmapped_minimap();
    }

    bool is_minimap_unmapped_floxer_mapped() const {
        return is_unmapped_minimap() && !is_unmapped_floxer();
    }

    bool is_floxer_unmapped_minimap_mapped() const {
        return !is_unmapped_minimap() && is_unmapped_floxer();
    }
};

void read_alignments(
    std::filesystem::path const& alignment_file_path,
    std::unordered_map<std::string, query_data_t>& query_data_by_query_id,
    bool const is_floxer,
    double const floxer_allowed_error_rate
) {
    seqan3::sam_file_input input{alignment_file_path};

    for (auto const& record : input) {
        auto iter = query_data_by_query_id.find(record.id());
        if (iter == query_data_by_query_id.end()) {
            auto emplace_result = query_data_by_query_id.emplace(record.id(), query_data_t{});
            assert(emplace_result.second);
            iter = emplace_result.first;
        }

        auto& query_data = iter->second;
        auto& alignment_data = is_floxer ? query_data.floxer_alignments : query_data.minimap_alignments;

        if (is_floxer) {
            query_data.mentioned_by_floxer = true;
        } else {
            query_data.mentioned_by_minimap = true;
        }

        if (static_cast<bool>(record.flag() & seqan3::sam_flag::unmapped)) {
            alignment_data.is_explicitly_unmapped = true;
        } else {
            alignment_data.is_mapped = true;
        }

        size_t num_query_bases_consumed_by_cigar = 0;
        size_t num_reference_bases_consumed_by_cigar = 0;
        size_t num_unclipped_query_bases_consumed_by_cigar = 0;
        size_t num_soft_clipped_bases = 0;
        size_t num_hard_clipped_bases = 0;

        size_t longest_indel = 0;

        for (auto const [count, operation] : record.cigar_sequence()) {
            if (
                operation == 'I'_cigar_operation || operation == 'M'_cigar_operation ||
                operation == '='_cigar_operation || operation == 'X'_cigar_operation
            ) {
                num_query_bases_consumed_by_cigar += count;
                num_unclipped_query_bases_consumed_by_cigar += count;
            }

            if (
                operation == 'D'_cigar_operation || operation == 'M'_cigar_operation ||
                operation == '='_cigar_operation || operation == 'X'_cigar_operation
            ) {
                num_reference_bases_consumed_by_cigar += count;
            }

            if (operation == 'S'_cigar_operation) {
                num_query_bases_consumed_by_cigar += count;
                num_soft_clipped_bases += count;
            }

            if (operation == 'H'_cigar_operation) {
                num_hard_clipped_bases += count;
            }

            if (operation == 'N'_cigar_operation) {
                spdlog::warn(
                    "Unexpected cigar operation in {} alignment of query {}: {}, count {}",
                    is_floxer ? "floxer" : "minimap",
                    record.id(),
                    operation.to_char(),
                    count
                );
            }

            if (operation == 'I'_cigar_operation || operation == 'D'_cigar_operation) {
                longest_indel = std::max(longest_indel, static_cast<size_t>(count));
            }
        }

        if (!record.sequence().empty() && num_hard_clipped_bases == 0) {
            if (query_data.sequence.has_value()) {
                if (
                    query_data.sequence.value() != record.sequence() &&
                    !std::ranges::equal(query_data.sequence.value(), record.sequence() | std::views::reverse | seqan3::views::complement)
                ) {
                    spdlog::warn(
                        "Observed different sequences for query {}. New one is by {}. Old length: {}, new length: {}",
                        record.id(),
                        is_floxer ? "floxer" : "minimap",
                        query_data.sequence.value().size(),
                        record.sequence().size()
                    );
                }
            } else {
                query_data.sequence = record.sequence();
            }
        }

        if (alignment_data.is_explicitly_unmapped) {
            continue;
        }

        size_t const query_length_without_clipped_bases = num_query_bases_consumed_by_cigar - num_hard_clipped_bases - num_soft_clipped_bases;
        size_t const edit_distance = static_cast<size_t>(record.tags().get<"NM"_tag>());
        double const edit_distance_error_rate = static_cast<double>(edit_distance) / static_cast<double>(query_length_without_clipped_bases);

        bool is_inversion = false;
        try {
            if (record.tags().get<"tp"_tag>() == 'I') {
                is_inversion = true;
            }
        } catch (std::out_of_range& exc) {
            // nothing to be done if the tag is not there
        }

        auto const extracted_record = alignment_record_t {
            .num_query_bases_consumed_by_cigar = num_query_bases_consumed_by_cigar,
            .num_unclipped_query_bases_consumed_by_cigar = num_unclipped_query_bases_consumed_by_cigar,
            .num_reference_bases_consumed_by_cigar = num_reference_bases_consumed_by_cigar,

            .num_soft_clipped_bases = num_soft_clipped_bases,
            .num_hard_clipped_bases = num_hard_clipped_bases,

            .query_length_without_clipped_bases = query_length_without_clipped_bases,

            .edit_distance = edit_distance,
            .edit_distance_error_rate = edit_distance_error_rate,

            .is_inversion = is_inversion,

            .longest_indel = longest_indel
        };

        if (
            !static_cast<bool>(record.flag() & seqan3::sam_flag::secondary_alignment) &&
            !static_cast<bool>(record.flag() & seqan3::sam_flag::supplementary_alignment)
        ) {
            if (alignment_data.primary_alignment.has_value()) {
                spdlog::warn("Multiple primary alignments for query");
            } else {
                alignment_data.primary_alignment.emplace(extracted_record);
            }

            continue;
        }

        if (static_cast<bool>(record.flag() & seqan3::sam_flag::supplementary_alignment)) {
            if (static_cast<bool>(record.flag() & seqan3::sam_flag::secondary_alignment)) {
                alignment_data.secondary_supplementary_alignments.emplace_back(extracted_record);
            }

            alignment_data.supplementary_alignments.emplace_back(extracted_record);
            continue;
        }

        if (extracted_record.is_inversion) {
            alignment_data.secondary_inverted_alignments.emplace_back(extracted_record);
        } else if (extracted_record.is_significantly_clipped(floxer_allowed_error_rate)) {
            alignment_data.secondary_linear_clipped_alignments.emplace_back(extracted_record);
        } else if (extracted_record.is_high_edit_distance(floxer_allowed_error_rate)) {
            alignment_data.secondary_linear_high_edit_distance_alignments.emplace_back(extracted_record);
        } else {
            alignment_data.secondary_linear_basic_alignments.emplace_back(extracted_record);
        }
    }
}

void validate_query_data(std::unordered_map<std::string, query_data_t> const& query_data_by_query_id) {
    for (auto const& [query_id, query_data] : query_data_by_query_id) {
        query_data.check_expectations(query_id);
    }
}

void print_value(
    std::string_view const predicate_name,
    size_t const num_queries_matching_predicate,
    [[maybe_unused]] double const num_subset_queries,
    [[maybe_unused]] double const num_queries
) {
    fmt::print(
        "{} = {}\n",
        predicate_name,
        num_queries_matching_predicate
    );
}

void print_basic_stats(std::unordered_map<std::string, query_data_t> const& query_data_by_query_id) {
    size_t const num_queries = query_data_by_query_id.size();

    size_t num_unmapped_floxer = 0;
    size_t num_unmapped_minimap = 0;
    size_t num_both_mapped = 0;
    size_t num_both_unmapped = 0;
    size_t num_minimap_unmapped_floxer_mapped = 0;
    size_t num_floxer_unmapped_minimap_mapped = 0;

    for (auto const& [query_id, query_data] : query_data_by_query_id) {
        if (query_data.is_unmapped_floxer()) {
            ++num_unmapped_floxer;
        }
        if (query_data.is_unmapped_minimap()) {
            ++num_unmapped_minimap;
        }

        if (query_data.is_mapped_both()) {
            ++num_both_mapped;
        }

        if (query_data.is_unmapped_both()) {
            ++num_both_unmapped;
        }

        if (query_data.is_floxer_unmapped_minimap_mapped()) {
            ++num_floxer_unmapped_minimap_mapped;
        }

        if (query_data.is_minimap_unmapped_floxer_mapped()) {
            ++num_minimap_unmapped_floxer_mapped;
        }
    }

    size_t const num_mapped_floxer = num_queries - num_unmapped_floxer;
    size_t const num_mapped_minimap = num_queries - num_unmapped_minimap;

    fmt::print("[general_stats]\n");
    print_value("number_of_queries", num_queries, num_queries, num_queries);
    print_value("both_mapped", num_both_mapped, num_queries, num_queries);
    print_value("both_unmapped", num_both_unmapped, num_queries, num_queries);
    print_value("floxer_mapped", num_mapped_floxer, num_queries, num_queries);
    print_value("floxer_unmapped", num_unmapped_floxer, num_queries, num_queries);
    print_value("minimap_mapped", num_mapped_minimap, num_queries, num_queries);
    print_value("minimap_unmapped", num_unmapped_minimap, num_queries, num_queries);
    print_value("floxer_unmapped_and_minimap_mapped", num_floxer_unmapped_minimap_mapped, num_queries, num_queries);
    print_value("minimap_unmapped_and_floxer_mapped", num_minimap_unmapped_floxer_mapped, num_queries, num_queries);
}

template<typename V>
concept alignment_data_view = std::ranges::view<V> &&
    std::same_as<std::ranges::range_value_t<V>, alignment_data_for_query_t>;

void print_alignment_statistics(
    std::string_view const title,
    size_t const num_queries,
    double const floxer_allowed_error_rate,
    alignment_data_view auto alignments
) {
    fmt::print("[{}]\n", title);

    size_t num_best_chimeric_or_inversion = 0;
    size_t num_best_significantly_clipped = 0;
    size_t num_best_high_edit_distance = 0;
    size_t num_basic = 0;

    size_t num_multiple_mapping = 0;

    double basic_longest_indel_sum = 0;
    double basic_alignments_error_rate_sum = 0.0;

    size_t num_subset_queries = 0;

    for (auto const& alignment_data : alignments) {
        if (alignment_data.best_is_chimeric_or_inversion(floxer_allowed_error_rate)) {
            ++num_best_chimeric_or_inversion;
        }

        if (alignment_data.best_is_significantly_clipped(floxer_allowed_error_rate)) {
            ++num_best_significantly_clipped;
        }

        if (alignment_data.best_is_high_edit_distance(floxer_allowed_error_rate)) {
            ++num_best_high_edit_distance;
        }

        if (alignment_data.has_basic(floxer_allowed_error_rate)) {
            basic_alignments_error_rate_sum += alignment_data.basic_error_rate(floxer_allowed_error_rate);
            basic_longest_indel_sum += alignment_data.basic_longest_indel(floxer_allowed_error_rate);
            ++num_basic;
        }

        if (alignment_data.is_multiple_mapping()) {
            ++num_multiple_mapping;
        }

        ++num_subset_queries;
    }

    // basic = linear, not significantly clipped, not higher error rate than in floxer config
    fmt::print("num_queries = {}\n", num_subset_queries);
    print_value("num_best_chimeric_or_inversion", num_best_chimeric_or_inversion, num_subset_queries, num_queries);
    print_value("num_best_significantly_clipped", num_best_significantly_clipped, num_subset_queries, num_queries);
    print_value("num_best_high_edit_distance", num_best_high_edit_distance, num_subset_queries, num_queries);
    print_value("num_basic", num_basic, num_subset_queries, num_queries);
    print_value("multiple_mapping", num_multiple_mapping, num_subset_queries, num_queries);
    fmt::print(
        "basic_average_longest_indel = {}\n",
        basic_longest_indel_sum / static_cast<double>(num_basic)
    );
    fmt::print(
        "basic_alignments_average_error_rate = {}\n",
        basic_alignments_error_rate_sum / static_cast<double>(num_basic)
    );
}

int main(int argc, char** argv) {
    sharg::parser parser{ "compare_aligner_outputs", argc, argv, sharg::update_notifications::off };

    parser.info.author = about_floxer::author;
    parser.info.description = {
        "Compares the alignment output of two readmappers regarding alignments found, edit distance and large indels. "
        "This program was created to compare specifically minimap2 and floxer."
    };
    parser.info.email = about_floxer::email;
    parser.info.url = about_floxer::url;
    parser.info.short_description = "Compare the alignment output of two readmappers.";
    parser.info.synopsis = {
        "./compare_aligner_outputs --reference minimap2_alignments.sam --new floxer_alignments.sam"
    };
    parser.info.version = "1.0.0";
    parser.info.date = "03.06.2024";

    std::filesystem::path minimap_input_path{};
    std::filesystem::path floxer_input_path{};
    double error_rate = 0.1;

    parser.add_option(minimap_input_path, sharg::config{
        .short_id = 'r',
        .long_id = "reference",
        .description = "The sam file of the reference read mapper (e.g. minimap2).",
        .required = true,
        .validator = sharg::input_file_validator{}
    });

    parser.add_option(floxer_input_path, sharg::config{
        .short_id = 'n',
        .long_id = "new",
        .description = "The sam file of the new read mapper (e.g. floxer).",
        .required = true,
        .validator = sharg::input_file_validator{}
    });

    parser.add_option(error_rate, sharg::config{
        .short_id = 'e',
        .long_id = "error-rate",
        .description = "The expected error rate of the aligners (especially floxer).",
        .required = false,
        .validator = sharg::arithmetic_range_validator{0.00001, 0.99999}
    });

    parser.parse();

    std::unordered_map<std::string, query_data_t> query_data_by_query_id{};

    read_alignments(minimap_input_path, query_data_by_query_id, false, error_rate);
    read_alignments(floxer_input_path, query_data_by_query_id, true, error_rate);

    validate_query_data(query_data_by_query_id);

    print_basic_stats(query_data_by_query_id);

    size_t const num_queries = query_data_by_query_id.size();

    print_alignment_statistics(
        "floxer_stats_if_floxer_mapped",
        num_queries,
        error_rate,
        query_data_by_query_id |
            std::views::values |
            std::views::filter([] (query_data_t const& query_data) {
                return !query_data.is_unmapped_floxer();
            }) |
            std::views::transform(&query_data_t::floxer_alignments)
    );
    print_alignment_statistics(
        "minimap_stats_if_minimap_mapped",
        num_queries,
        error_rate,
        query_data_by_query_id |
            std::views::values |
            std::views::filter([] (query_data_t const& query_data) {
                return !query_data.is_unmapped_minimap();
            }) |
            std::views::transform(&query_data_t::minimap_alignments)
    );
    print_alignment_statistics(
        "minimap_stats_if_both_mapped",
        num_queries,
        error_rate,
        query_data_by_query_id |
            std::views::values |
            std::views::filter([] (query_data_t const& query_data) {
                return query_data.is_mapped_both();
            }) |
            std::views::transform(&query_data_t::minimap_alignments)
    );
    print_alignment_statistics(
        "minimap_stats_if_only_minimap_mapped",
        num_queries,
        error_rate,
        query_data_by_query_id |
            std::views::values |
            std::views::filter([] (query_data_t const& query_data) {
                return query_data.is_floxer_unmapped_minimap_mapped();
            }) |
            std::views::transform(&query_data_t::minimap_alignments)
    );

    return 0;
}
