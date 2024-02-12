#include <cli.hpp>
#include <io.hpp>
#include <pex.hpp>
#include <search.hpp>

#include <algorithm>
#include <unordered_map>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <fmindex-collection/BiFMIndex.h>
#include <fmindex-collection/occtable/InterleavedEPRV2.h>
#include <fmindex-collection/search/SearchNg21.h>

int main(int argc, char** argv) {
    auto const opt = cli::parse_options(argc, argv);

    fmt::println("/\\ \\/ /\\ \\/ /\\ welcome to floxer /\\ \\/ /\\ \\/ /\\");

    auto const input = io::read_inputs(opt.reference_sequence, opt.queries);

    size_t constexpr Sigma = 5; // DNA + Sentinel
    using Table = fmindex_collection::occtable::interleavedEPR16V2::OccTable<Sigma>;

    size_t const sampling_rate = 16; // what are good values for my use case?
    size_t const num_threads = 1;

    auto const index = fmindex_collection::BiFMIndex<Table>(
        input.reference_sequences, sampling_rate, num_threads
    );

    // for now assume every leaf has the same number of errors
    search_scheme_cache scheme_cache(opt.pex_leaf_num_errors);
    pex_tree_cache tree_cache{};

    for (auto const& query : input.queries) {
        auto const tree_config = pex_tree_config {
            .total_query_length = query.sequence.size(),
            .query_num_errors = opt.query_num_errors,
            .leaf_num_errors = opt.pex_leaf_num_errors
        };

        auto const& pex_tree = tree_cache.get(tree_config);
        // for now assume every leaf has the same length
        auto const& search_scheme = scheme_cache.get(pex_tree.leaf_query_length());
        auto const leaf_queries = pex_tree.generate_leaf_queries(query.sequence);

        std::vector<std::unordered_map<size_t, size_t>> useful_hits(leaf_queries.size());

        fmindex_collection::search_ng21::search(
            index, leaf_queries, search_scheme, [&index, &useful_hits] (size_t const query_id, auto cursor, size_t const errors) {                
                auto & useful_query_hits = useful_hits[query_id];

                for (auto hits{begin(cursor)}; hits < end(cursor); ++hits) {
                    auto const [reference_id, pos] = index.locate(hits);

                    size_t const start = errors <= pos ? pos - errors : 0ul; // protect against underflow
                    size_t const end = pos + errors;
                    bool better_hit_exists = false;
                    for (size_t neighborhood_pos = start; neighborhood_pos <= end; ++neighborhood_pos) {
                        auto const iter = useful_query_hits.find(neighborhood_pos);
                        if (iter != useful_query_hits.end() && iter->second < errors) {
                            better_hit_exists = true;
                            break;
                        }
                    }

                    // if we can't assume that the hits are reported in asceding order of number of errors,
                    // some cleanup needs to be done here or later

                    if (better_hit_exists) {
                        continue;
                    }

                    useful_query_hits[pos] = errors;
                }
            }
        );

        for (size_t i = 0; i < leaf_queries.size(); ++i) {
            fmt::println("query: {}", leaf_queries[i]);
            for (auto const & entry : useful_hits[i]) {
                fmt::println("    - {}, {} error{}", entry.first, entry.second, entry.second == 1 ? "" : "s");
            }
        }
    }

    return 0;
}
