#include <gtest/gtest.h>

#include <fmindex.hpp>
#include <pex.hpp>
#include <search.hpp>

TEST(floxer_test, BasicAssertions) {
    std::vector<std::vector<uint8_t>> references{{
        0,0,0,0,0,0,0,0,0,0,
        1,1,1,1,1,1,1,1,1,1,
        2,2,2,2,2,2,2,2,2,2,
        3,3,3,3,3,3,3,3,3,3
    }};

    fmindex index(references, 16, 1);
    
    // TODO finish this test once the search interface is clear.
    // Right now it would be a copy and paste of the floxxer main function.
}
