#include <catch2/catch_all.hpp>
#include <iostream>
#include "utils/combinations.hpp"
#include "utils/utils.hpp"

using std::vector;

TEST_CASE("combinations-k=2", "[utils]") {

    std::vector<int> c({1, 2, 3, 4});
    auto combs = ggbp::utils::Combinations(c, 2);

    std::vector<std::vector<int>> expected = {
        {1, 2},
        {1, 3},
        {1, 4},
        {2, 3},
        {2, 4},
        {3, 4}
    };

    int i = 0;
    for (auto s : combs) {
        if (i >= 6) REQUIRE(false);
        REQUIRE(s == expected[i]);
        ++i;
    }
}


TEST_CASE("combinations-k=1", "[utils]") {

    std::vector<int> c({1, 2, 3, 4});
    auto combs = ggbp::utils::Combinations(c, 1);

    std::vector<std::vector<int>> expected = {
        {1},
        {2},
        {3},
        {4},
    };

    int i = 0;
    for (auto s : combs) {
        if (i >= 4) REQUIRE(false);
        REQUIRE(s == expected[i]);
        ++i;
    }
}


TEST_CASE("combinations-k=0", "[utils]") {

    std::vector<int> c({1, 2, 3, 4});
    auto combs = ggbp::utils::Combinations(c, 0);

    int i = 0;
    for (auto s : combs) {
        if (i >= 1) REQUIRE(false);
        REQUIRE(s.size() == 0);
        ++i;
    }
}

TEST_CASE("combinations-k=4", "[utils]") {

    std::vector<int> c({1, 2, 3, 4});
    auto combs = ggbp::utils::Combinations(c, 4);

    std::vector<std::vector<int>> expected = {
        {1, 2, 3, 4},
    };

    int i = 0;
    for (auto s : combs) {
        if (i >= 1) REQUIRE(false);
        REQUIRE(s == expected[i]);
        ++i;
    }
}
