#include <cli.hpp>
#include <io.hpp>
#include <pex.hpp>

#include <iostream>

// temporary
#include <span>

#include <fmindex-collection/BiFMIndex.h>
#include <fmindex-collection/occtable/InterleavedEPRV2.h>
#include <fmindex-collection/search/SearchNg21.h>
#include <search_schemes/generator/optimum.h>

int main(int argc, char** argv) {
    auto const opt = cli::parse_options(argc, argv);

    std::cout << "---{ welcome to floxer }---\n\n"
        << "- reference path: " << opt.reference_genome.c_str() << '\n'
        << "- query path: " << opt.queries.c_str() << '\n'
        << "- output path: " << opt.output_file.c_str() << '\n'
        << "- number of allowed errors in query: " << opt.query_num_errors << '\n'
        << "- number of errors in PEX leaves: " << opt.pex_leaf_num_errors << "\n\n";

    auto const input = io::read_inputs(opt.reference_genome, opt.queries);

    size_t constexpr Sigma = 5; // DNA + Sentinel
    using Table = fmindex_collection::occtable::interleavedEPR16V2::OccTable<Sigma>;

    size_t const sampling_rate = 16; // what are good values for my use case?
    size_t const num_threads = 1;

    std::vector<std::vector<uint8_t>> const reference_genomes{std::move(input.reference_genome)};
    auto const index = fmindex_collection::BiFMIndex<Table>(
        reference_genomes, sampling_rate, num_threads
    );

    // SIMON is the optimum search scheme the one I want?
    auto const search_scheme = search_schemes::generator::optimum(0, opt.query_num_errors);

    pex_tree_cache tree_cache{};
    for (auto const& query : input.queries) {
        auto const tree_config = pex_tree_config {
            .total_query_length = query.sequence.size(),
            .query_num_errors = opt.query_num_errors,
            .leaf_num_errors = opt.pex_leaf_num_errors
        };

        auto const& pex_tree = tree_cache.get(tree_config);

        std::vector<std::span<const uint8_t>> leaf_queries(pex_tree.get_leaves().size());
        for (auto const& leaf : pex_tree.get_leaves()) {
            size_t const length = leaf.query_index_to - leaf.query_index_from + 1;
            auto const leaf_query_span = std::span(query.sequence).subspan(leaf.query_index_from, length);
            // workaround copy - why does is break down - this caused a segfault
            leaf_queries.emplace_back(std::move(leaf_query_span)/*.begin(), leaf_query_span.end()*/);
        }   

        std::cout << "Searching leaves for query " << query.tag << std::endl;

        fmindex_collection::search_ng21::search(
            index, leaf_queries, search_scheme, [&index] (size_t query_id, auto cursor, size_t errors) {
                std::cout << "found leaf " << query_id 
                    << ", " << cursor.count() << " times and " 
                    << errors << " errors" << std::endl; 

                for (auto i{begin(cursor)}; i < end(cursor); ++i) {
                    auto const [reference_id, pos] = index.locate(i);
                    std::cout << "-reference " << reference_id << " position " << pos << std::endl;
                }
            }
        );
    }

    std::cout << "---{ Verifiying hits..." << std::endl;

    // TODO

    std::cout << "                              ...done }---\n\n";

    return 0;
}
