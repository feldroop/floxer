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

bool better_hit_exists_at(
    size_t const pos, 
    size_t const threshold,
    std::unordered_map<size_t, size_t> & hits
) {
    auto const iter = hits.find(pos);
    return iter != hits.end() && iter->second <= threshold;
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
    auto const queries = io::read_queries(opt.queries);

    // FIX LATER for now assume every leaf has the same number of errors
    search_scheme_cache scheme_cache(opt.pex_leaf_num_errors);
    pex_tree_cache tree_cache{};

    for (auto const& query : queries) {
        fmt::println("query {}:", query.tag);

        auto const tree_config = pex_tree_config {
            .total_query_length = query.sequence.size(),
            .query_num_errors = opt.query_num_errors,
            .leaf_num_errors = opt.pex_leaf_num_errors
        };

        auto const& pex_tree = tree_cache.get(tree_config);
        // FIX LATER for now assume every leaf has the same length
        auto const& search_scheme = scheme_cache.get(pex_tree.leaf_query_length());
        auto const leaf_queries = pex_tree.generate_leaf_queries(query.sequence);
        
        // REFACTOR LATER useful_hits[query_id][reference_id][hit_pos] -> num_errors
        using hit_map = std::vector<std::vector<std::unordered_map<size_t, size_t>>>;
        hit_map useful_hits(
            leaf_queries.size(), 
            std::vector<std::unordered_map<size_t, size_t>>(num_reference_sequences)
        );

        fmindex_collection::search_ng21::search(
            index,
            leaf_queries,
            search_scheme,
            [&index, &useful_hits] (size_t const query_id, auto cursor, size_t const errors) {                
                auto & useful_query_hits = useful_hits[query_id];

                for (auto hit{begin(cursor)}; hit < end(cursor); ++hit) {
                    auto const [reference_id, pos] = index.locate(hit);
                    auto & useful_query_to_reference_hits = useful_query_hits[reference_id];

                    for (size_t dist = errors; dist != std::numeric_limits<size_t>::max(); --dist) {
                        size_t const threshold = errors - dist;
                        size_t const upper_pos = pos + dist;
                        if (better_hit_exists_at(upper_pos, threshold, useful_query_to_reference_hits)) {
                            break;
                        }

                        // this if is only true if we checked all of the possible better hits
                        // hence we insert now and break to not do the repetitive check
                        if (dist == 0) {
                            // FIX LATER if we can't assume that the hits are reported in asceding order of 
                            // number of errors, some cleanup needs to be done here or later
                            useful_query_to_reference_hits[pos] = errors;
                            break;
                        }

                        size_t const lower_pos = dist <= pos ? pos - dist : 0ul; // protect against underflow
                        if (better_hit_exists_at(lower_pos, threshold, useful_query_to_reference_hits)) {
                            break;
                        }
                    }
                }
            }
        );

        for (size_t query_id = 0; query_id < leaf_queries.size(); ++query_id) {
            fmt::println("    leaf: {}", leaf_queries[query_id]);
            for (size_t reference_id = 0; reference_id < num_reference_sequences; ++ reference_id) {
                for (auto const & entry : useful_hits[query_id][reference_id]) {
                    fmt::println(
                        "        - at {}, reference: {}, {} error{}", 
                        entry.first,
                        index_and_data.reference_tags[reference_id],
                        entry.second,
                        entry.second == 1 ? "" : "s"
                    );
                }
            }
        }
    }

    return 0;
}
