#include <cli.hpp>
#include <io.hpp>
#include <pex.hpp>
#include <search.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <fmindex-collection/BiFMIndex.h>
#include <fmindex-collection/occtable/InterleavedEPRV2.h>
#include <fmindex-collection/search/SearchNg22.h>

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

    search_scheme_cache scheme_cache(opt.pex_leaf_num_errors);
    pex_tree_cache tree_cache{};

    for (auto const& query : input.queries) {
        auto const tree_config = pex_tree_config {
            .total_query_length = query.sequence.size(),
            .query_num_errors = opt.query_num_errors,
            .leaf_num_errors = opt.pex_leaf_num_errors
        };

        auto const& pex_tree = tree_cache.get(tree_config);
        auto const& search_scheme = scheme_cache.get(pex_tree.leaf_query_length());
        auto const leaf_queries = pex_tree.generate_leaf_queries(query.sequence);

        fmt::print(
            "queries:\n{}\nquery length: {}\n", 
            leaf_queries, 
            pex_tree.leaf_query_length()
        );

        fmindex_collection::search_ng22::search(
            index, leaf_queries, search_scheme, [&index] (size_t const query_id, auto cursor, size_t const errors, auto const& thing) {
                fmt::println("found leaf {}, {} times, {} errors, alignment: {}", query_id, cursor.count(), errors, thing);
                
                for (auto i{begin(cursor)}; i < end(cursor); ++i) {
                    auto const [reference_id, pos] = index.locate(i);
                    fmt::println("- reference {}, position {}", reference_id, pos);
                }
            }
        );
    }

    return 0;
}
