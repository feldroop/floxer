#include <cli.hpp>
#include <fmindex.hpp>
#include <input.hpp>
#include <output.hpp>
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
    
    cli::options opt;
    try {
        opt = cli::parse_and_validate_options(argc, argv);
    } catch (std::exception const & e) {
        fmt::print(stderr, "[CLI PARSER ERROR]\n{}\n", e.what());
        return -1;
    }
    
    std::vector<input::reference_record> references;
    try {
        references = input::read_references(opt.reference_sequence);
    } catch (std::exception const& e) {
        fmt::print(
            stderr,
            "[INPUT ERROR]\nAn error occured while trying to read the reference from "
            "the file {}.\n{}\n",
            opt.reference_sequence.c_str(),
            e.what()
        );
        return -1;
    }

    fmindex index;
    if (!opt.index_path.empty() && std::filesystem::exists(opt.index_path)) {
        try {
            index = input::load_index(opt.index_path);
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
        size_t const suffix_array_sampling_rate = 16; 

        index = fmindex(
            references | std::views::transform(&input::reference_record::sequence),
            suffix_array_sampling_rate,
            opt.num_threads
        );

        if (!opt.index_path.empty()) {
            try {
                output::save_index(index, opt.index_path);
            } catch (std::exception const& e) {
                fmt::print(
                    stderr,
                    "[OUTPUT WARNING]\nAn error occured while trying to write the index to "
                    "the file {}.\nContinueing without saving the index.\n{}\n",
                    opt.index_path.c_str(),
                    e.what()
                );
            }
        }
    }

    std::vector<input::query_record> fastq_queries;
    try {
        fastq_queries = input::read_queries(opt.queries);
    } catch (std::exception const& e) {
        fmt::print(
            stderr,
            "[INPUT ERROR]\nAn error occured while trying to read the queries from "
            "the file {}.\n{}\n",
            opt.queries.c_str(),
            e.what()
        );
        return -1;
    }

    auto const output_header = output::sam_header(references);

    search::search_scheme_cache scheme_cache;
    pex_tree_cache tree_cache;

    for (auto const& fastq_query : fastq_queries) {
        size_t const query_num_errors = fastq_query.num_errors_from_user_config(opt);

        auto const tree_config = pex_tree_config {
            .total_query_length = fastq_query.sequence_length,
            .query_num_errors = query_num_errors,
            .leaf_max_num_errors = opt.pex_leaf_num_errors
        };

        auto const& tree = tree_cache.get(tree_config);
        auto const alignments = tree.search(
            references,
            fastq_query.sequence,
            scheme_cache,
            index
        );

        bool found_any_alignments = false;
        // TODO create and write .sam record here
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
