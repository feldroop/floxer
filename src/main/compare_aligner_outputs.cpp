#include <about_floxer.hpp>

#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <sharg/all.hpp>
#include <spdlog/spdlog.h>

using namespace seqan3::literals;

// this needs to be added, because the tp tag is used defined by minimap
template <>
struct seqan3::sam_tag_type<"tp"_tag> {
    using type = char;
};

struct alignment_data_t {
    bool is_unmapped_floxer{true};
    bool is_unmapped_minimap{true};

    bool is_non_linear_minimap{false};

    bool mentioned_by_floxer{false};
    bool mentioned_by_minimap{false};

    size_t minimap_longest_indel{0};
    size_t floxer_longest_indel{0};
};

void read_alignments(
    std::filesystem::path const& alignment_file_path,
    std::unordered_map<std::string, alignment_data_t>& alignment_data_by_query_id,
    bool const is_floxer,
    double const error_rate
) {
    seqan3::sam_file_input input{alignment_file_path};

    for (auto const& record : input) {
        auto iter = alignment_data_by_query_id.find(record.id());
        if (iter == alignment_data_by_query_id.end()) {
            auto emplace_result = alignment_data_by_query_id.emplace(record.id(), alignment_data_t{});
            assert(emplace_result.second);
            iter = emplace_result.first;
        }

        auto& alignment_data = iter->second;

        if (is_floxer) {
            alignment_data.mentioned_by_floxer = true;
        } else {
            alignment_data.mentioned_by_minimap = true;
        }

        if (static_cast<bool>(record.flag() & seqan3::sam_flag::unmapped)) {
            continue;
        }

        if (static_cast<bool>(record.flag() & seqan3::sam_flag::supplementary_alignment)) {
            if (is_floxer) {
                spdlog::warn("Unexpected non-linear floxer alignment");
            }
            // is chimeric
            alignment_data.is_non_linear_minimap = true;
        }

        if (!is_floxer) {
            try {
                if (record.tags().get<"tp"_tag>() == 'I') {
                    // is inversion
                    alignment_data.is_non_linear_minimap = true;
                }
            } catch (std::out_of_range& exc) {
                // nothing to be done if the tag is not there
            }
        }

        size_t cigar_length = 0;
        for (auto const [count, operation] : record.cigar_sequence()) {
            cigar_length += count;

            if (operation == 'I'_cigar_operation || operation == 'D'_cigar_operation) {
                if (is_floxer) {
                    alignment_data.floxer_longest_indel = std::max(alignment_data.floxer_longest_indel, static_cast<size_t>(count));
                } else {
                    alignment_data.minimap_longest_indel = std::max(alignment_data.minimap_longest_indel, static_cast<size_t>(count));
                }
            }
        }

        size_t const max_num_errors = cigar_length * error_rate + 1;
        if (record.tags().get<"NM"_tag>() > static_cast<int>(max_num_errors)) {
            if (is_floxer) {
                spdlog::warn(
                    "Unexpected floxer alignment with large number of errors. Size {}, expected max errors: {}, actual: {}",
                    cigar_length,
                    max_num_errors,
                    record.tags().get<"NM"_tag>()
                );
            } else {
                alignment_data.is_non_linear_minimap = true;
            }
        }

        if (is_floxer) {
            alignment_data.is_unmapped_floxer = false;
        } else {
            alignment_data.is_unmapped_minimap = false;
        }
    }
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

    std::unordered_map<std::string, alignment_data_t> alignment_data_by_query_id{};

    read_alignments(minimap_input_path, alignment_data_by_query_id, false, error_rate);
    read_alignments(floxer_input_path, alignment_data_by_query_id, true, error_rate);

    size_t const num_queries = alignment_data_by_query_id.size();

    size_t num_unmapped_floxer = 0;
    size_t num_unmapped_minimap = 0;

    size_t num_non_linear_minimap = 0;

    size_t num_unmapped_both = 0;
    size_t num_minimap_unmapped_floxer_mapped = 0;
    size_t num_floxer_unmapped_minimap_linear_mapped = 0;
    size_t num_floxer_unmapped_minimap_non_linear_mapped = 0;
    size_t num_mapped_both_minimap_linear = 0;
    size_t num_mapped_both_minimap_non_linear = 0;

    size_t minimap_longest_indel_sum = 0;
    size_t floxer_longest_indel_sum = 0;

    size_t floxer_unmapped_minimap_mapped_minimap_longest_indel_sum = 0;

    for (auto const& [query_id, alignment_data] : alignment_data_by_query_id) {
        if (!alignment_data.mentioned_by_floxer) {
            spdlog::warn("Query {} not mentioned by floxer", query_id);
        }

        if (!alignment_data.mentioned_by_minimap) {
            spdlog::warn("Query {} not mentioned by minimap", query_id);
        }

        if (alignment_data.is_unmapped_floxer) {
            ++num_unmapped_floxer;
        }

        if (alignment_data.is_unmapped_minimap) {
            ++num_unmapped_minimap;
        }

        if (alignment_data.is_non_linear_minimap) {
            ++num_non_linear_minimap;
        }

        if (alignment_data.is_unmapped_floxer && alignment_data.is_unmapped_minimap) {
            ++num_unmapped_both;
        }

        if (!alignment_data.is_unmapped_floxer && alignment_data.is_unmapped_minimap) {
            ++num_minimap_unmapped_floxer_mapped;
        }

        if (alignment_data.is_unmapped_floxer && !alignment_data.is_unmapped_minimap) {
            floxer_unmapped_minimap_mapped_minimap_longest_indel_sum += alignment_data.minimap_longest_indel;

            if (alignment_data.is_non_linear_minimap) {
                ++num_floxer_unmapped_minimap_non_linear_mapped;
            } else {
                ++num_floxer_unmapped_minimap_linear_mapped;
            }
        }

        if (!alignment_data.is_unmapped_floxer && !alignment_data.is_unmapped_minimap) {
            if (alignment_data.is_non_linear_minimap) {
                ++num_mapped_both_minimap_non_linear;
            } else {
                ++num_mapped_both_minimap_linear;
            }
        }

        minimap_longest_indel_sum += alignment_data.minimap_longest_indel;
        floxer_longest_indel_sum  += alignment_data.floxer_longest_indel;
    }

    double const num_queries_d = static_cast<double>(num_queries);
    double const num_floxer_unmapped_minimap_mapped_d = static_cast<double>(
        num_floxer_unmapped_minimap_linear_mapped + num_floxer_unmapped_minimap_non_linear_mapped
    );

    spdlog::info("Number of queries: {} ({:.2f})", num_queries, num_queries / num_queries_d);
    spdlog::info("Floxer unmapped queries: {} ({:.2f})", num_unmapped_floxer, num_unmapped_floxer / num_queries_d);
    spdlog::info("Minimap unmapped queries: {} ({:.2f})", num_unmapped_minimap, num_unmapped_minimap / num_queries_d);
    spdlog::info(
        "Minimap non-linear mapped queries: {} ({:.2f})",
        num_non_linear_minimap,
        num_non_linear_minimap / num_queries_d
    );
    spdlog::info("Both unmapped: {} ({:.2f})", num_unmapped_both, num_unmapped_both / num_queries_d);
    spdlog::info(
        "Floxer mapped, minimap unmapped: {} ({:.2f})",
        num_minimap_unmapped_floxer_mapped,
        num_minimap_unmapped_floxer_mapped / num_queries_d
    );
    spdlog::info(
        "Floxer unmapped, minimap linear mapped: {} ({:.2f})",
        num_floxer_unmapped_minimap_linear_mapped,
        num_floxer_unmapped_minimap_linear_mapped / num_queries_d
    );
    spdlog::info(
        "Floxer unmapped, minimap non-linear mapped: {} ({:.2f})",
        num_floxer_unmapped_minimap_non_linear_mapped,
        num_floxer_unmapped_minimap_non_linear_mapped / num_queries_d
    );
    spdlog::info(
        "Both mapped, minimap linear: {} ({:.2f})",
        num_mapped_both_minimap_linear,
        num_mapped_both_minimap_linear / num_queries_d
    );
    spdlog::info(
        "Both mapped, minimap non-linear: {} ({:.2f})",
        num_mapped_both_minimap_non_linear,
        num_mapped_both_minimap_non_linear / num_queries_d
    );
    spdlog::info("floxer average longest indel: {:.2f}", floxer_longest_indel_sum / num_queries_d);
    spdlog::info("minimap average longest indel: {:.2f}", minimap_longest_indel_sum / num_queries_d);
    spdlog::info(
        "minimap average longest indel (floxer unmapped, minimap mapped): {:.2f}",
        floxer_unmapped_minimap_mapped_minimap_longest_indel_sum / num_floxer_unmapped_minimap_mapped_d
    );

    return 0;
}
