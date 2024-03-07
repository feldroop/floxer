#include <cli.hpp>
#include <fmindex.hpp>
#include <io.hpp>
#include <pex.hpp>
#include <search.hpp>

#include <exception>
#include <filesystem>
#include <ranges>
#include <span>
#include <vector>

#include <fmt/core.h>
#include <fmt/ranges.h>

int main(int argc, char** argv) {
    fmt::println("/\\ \\/ /\\ \\/ /\\ welcome to floxer /\\ \\/ /\\ \\/ /\\");

    auto const opt = cli::parse_and_validate_options(argc, argv);

    auto const references = io::read_references(opt.reference_sequence);
    fmindex index;
    if (!opt.index_path.empty() && std::filesystem::exists(opt.index_path)) {
        try {
            index = io::load_index(opt.index_path);
        } catch (std::exception const& e) {
            fmt::print(
                stderr,
                "[INPUT ERROR]\nAn error occured while trying to load the index from "
                "the file {}.\n{}\n",
                opt.index_path.c_str(),
                e.what()
            );
            return -1;
        }
    } else {
        // FIGURE OUT LATER what are good values for my use case?
        size_t const suffix_array_sampling_rate = 16; 

        index = fmindex(
            references | std::views::transform(&io::record::sequence),
            suffix_array_sampling_rate,
            opt.num_threads
        );

        if (!opt.index_path.empty()) {
            io::save_index(index, opt.index_path);
        }
    }

    auto const fastq_queries = io::read_queries(opt.queries);

    // FIX LATER for now assume every leaf has the same number of errors
    search::search_scheme_cache scheme_cache(opt.pex_leaf_num_errors);
    pex_tree_cache tree_cache{};

    for (auto const& fastq_query : fastq_queries) {
        fmt::println("query {}:", fastq_query.tag);

        auto const tree_config = pex_tree_config {
            .total_query_length = fastq_query.sequence.size(),
            .query_num_errors = opt.query_num_errors,
            .leaf_num_errors = opt.pex_leaf_num_errors
        };

        auto const& tree = tree_cache.get(tree_config);
        auto const alignments = tree.search(
            references,
            fastq_query.sequence,
            scheme_cache,
            index
        );

        bool found_any_alignments = false;
        for (size_t reference_id = 0; reference_id < alignments.size(); ++reference_id) {
            auto const& reference_alignments = alignments[reference_id];
            
            if (!reference_alignments.empty()) {
                fmt::println("\tto reference {}", reference_id);
                found_any_alignments = true;
            }

            for (auto const& alignment : std::views::values(reference_alignments)) {
                fmt::println("\t\t- {}", alignment);
            }
        }

        if (!found_any_alignments) {
            fmt::println("\t found no alignments");
        }
    }

    return 0;
}
