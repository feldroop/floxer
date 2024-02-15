#include <cli.hpp>
#include <floxer_fmindex.hpp>
#include <io.hpp>
#include <pex.hpp>
#include <search.hpp>

#include <algorithm>
#include <limits>
#include <unordered_map>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <fmindex-collection/search/SearchNg21.h>

struct hit {
    size_t position;
    size_t num_errors;

    bool is_better_than(hit const& other) {
        size_t const position_difference = position < other.position ?
            other.position - position : position - other.position;

        return num_errors <= other.num_errors && 
            position_difference <= other.num_errors - num_errors;
    }
};

// REFACTOR LATER hits[leaf_query_id][reference_id] -> hits
using hit_list = std::vector<std::vector<std::vector<hit>>>;

void erase_useless_hits(hit_list & hits) {
    static constexpr size_t erase_marker = std::numeric_limits<size_t>::max();
    
    for (auto & query_hits : hits) {
        for (auto & query_to_ref_hits : query_hits) {
            // this must stay, otherwise the expression in the below for loop head could underflow
            if (query_to_ref_hits.empty()) {
                continue;
            }

            std::ranges::sort(query_to_ref_hits, {}, [] (hit const& h) { return h.position; });

            for (size_t i = 0; i < query_to_ref_hits.size() - 1; ++i) {
                auto & current_hit = query_to_ref_hits[i];
                size_t j = i + 1;

                while (j < query_to_ref_hits.size() && current_hit.is_better_than(query_to_ref_hits[j])) { 
                    query_to_ref_hits[j].num_errors = erase_marker;
                    ++j; 
                }

                if (j < query_to_ref_hits.size() && query_to_ref_hits[j].is_better_than(current_hit)) {
                    current_hit.num_errors = erase_marker;
                }
            }

            std::erase_if(query_to_ref_hits, [] (hit const& h) { return h.num_errors == erase_marker; } );
        }
    }
}

hit_list search_fastq_query(
    io::query const& fastq_query,
    FloxerFMIndex const& index,
    pex_tree const& tree,
    search_scheme_cache& scheme_cache,
    size_t const num_reference_sequences
) {
    // FIX LATER for now assume every leaf has the same length
    auto const& search_scheme = scheme_cache.get(tree.leaf_query_length());
    auto const leaf_queries = tree.generate_leaf_queries(fastq_query.sequence);
    
    hit_list hits(
        leaf_queries.size(), 
        std::vector<std::vector<hit>>(num_reference_sequences)
    );

    fmindex_collection::search_ng21::search(
        index,
        leaf_queries,
        search_scheme,
        [&index, &hits] (size_t const leaf_query_id, auto cursor, size_t const errors) {                
            auto& query_hits = hits[leaf_query_id];

            for (auto hit{begin(cursor)}; hit < end(cursor); ++hit) {
                auto const [reference_id, pos] = index.locate(hit);
                query_hits[reference_id].emplace_back(pos, errors);
            }
        }
    );

    erase_useless_hits(hits);

    return hits;
}

void print_hits(
    hit_list const& hits,
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

    FloxerFMIndexWithMetaData index_and_data;
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
        }
    } else {
        auto reference_input = io::read_reference(opt.reference_sequence);
        
        size_t const suffix_array_sampling_rate = 16; // FIGURE OUT LATER what are good values for my use case?
        
        index_and_data.index = FloxerFMIndex(
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
    search_scheme_cache scheme_cache(opt.pex_leaf_num_errors);
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
            fastq_query,
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
