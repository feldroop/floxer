#include <cli.hpp>
#include <fmindex.hpp>
#include <io.hpp>
#include <pex.hpp>
#include <search.hpp>

#include <exception>
#include <vector>
#include <span>

#include <fmt/core.h>
#include <fmt/ranges.h>

void print_hits(
    search::hit_list const& hits,
    std::vector<std::span<const uint8_t>> leaf_queries,
    std::vector<std::string> const& reference_tags
) {
    for (size_t leaf_query_id = 0; leaf_query_id < hits.size(); ++leaf_query_id) {
        fmt::println("    leaf: {}", leaf_queries[leaf_query_id]);

        auto const& query_hits = hits[leaf_query_id];
        for (size_t reference_id = 0; reference_id < query_hits.size(); ++ reference_id) {
            for (auto const& hit : query_hits[reference_id]) {
                fmt::println(
                    "        - at {}, reference: {}, {} error{}", 
                    hit.position,
                    reference_tags[reference_id],
                    hit.num_errors,
                    hit.num_errors == 1 ? "" : "s"
                );
            }
        }
    }
}

int main(int argc, char** argv) {
    fmt::println("/\\ \\/ /\\ \\/ /\\ welcome to floxer /\\ \\/ /\\ \\/ /\\");

    auto const opt = cli::parse_and_validate_options(argc, argv);

    fmindex_with_metadata index_and_data;
    if (opt.reference_sequence.empty()) {
        try {
            index_and_data = io::load_index_and_data(opt.index_path);
        } catch (const std::exception& e) {
            fmt::print(
                stderr,
                "[INPUT ERROR]\nAn error occured while trying to load the index from "
                "the file {}.\n{}\n",
                opt.index_path.c_str(),
                e.what()
            );
            exit(-1);
        }
    } else {
        auto reference_input = io::read_reference(opt.reference_sequence);

        // FIGURE OUT LATER what are good values for my use case?
        size_t const suffix_array_sampling_rate = 16; 

        index_and_data.index = fmindex(
            reference_input.sequences,
            suffix_array_sampling_rate,
            opt.num_threads
        );

        index_and_data.reference_tags = std::move(reference_input.tags);

        if (!opt.index_path.empty()) {
            io::save_index_and_data(index_and_data, opt.index_path);
        }
    }

    size_t const num_reference_sequences = index_and_data.reference_tags.size();
    auto const& index = index_and_data.index;
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

        auto const hits = search_fastq_query(
            fastq_query.sequence,
            index,
            tree, 
            scheme_cache,
            num_reference_sequences
        );

        // temporary
        print_hits(
            hits, 
            tree.generate_leaf_queries(fastq_query.sequence), 
            index_and_data.reference_tags
        );
    }

    return 0;
}
