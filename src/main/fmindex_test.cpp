#include <fmindex.hpp>
#include <vector>

#include <fmindex-collection/search/SearchNg22.h>
#include <search_schemes/generator/optimum.h>
#include <search_schemes/expand.h>

#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/ranges.h>

int main() {
    std::vector<std::vector<uint8_t>> t{{2,2,2,1,1,1,1,2,2,2}};

    fmindex index(
        t,
        1,
        1
    );

    std::vector<std::vector<uint8_t>> q {{1,1,1,1}};

    auto search_scheme = search_schemes::expand(
            search_schemes::generator::optimum(0, 2),
        q[0].size()
    );

    fmindex_collection::search_ng22::search(
        index,
        q,
        search_scheme,
        [&index] (
            [[maybe_unused]] size_t const query_id,
            auto cursor,
            size_t const errors,
            auto alignment
        ) {
            fmt::print("{} errors, {}\n", errors, alignment);

            for (auto const& anchor : cursor) {
                auto const [reference_id, position] = index.locate(anchor);
                fmt::print("\tpos: {}\n", position);
            }
        }
    );

    fmt::print("done\n");
}
