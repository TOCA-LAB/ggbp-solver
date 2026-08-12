#include <catch2/catch_all.hpp>
#include <iostream>
#include <vector>

#include "ggbp/graph.hpp"
#include "utils/fix_c++_print.hpp"

#include "ggbp/ggbp.hpp"
#include "utils/utils.hpp"

TEST_CASE("get_best_ignitor", "[greedy]") {
    
    // Path with 10 vertices
    Graph G("IhCGGC@?G");

    std::vector<int> expected = {1, 5, 9, 2};

    std::vector<Vertex> burning_seq = ggbp::heuristics::spontaneous_neighbors(G, 1);

    for (auto [i, v] : ggbp::utils::enumerate(expected))
      REQUIRE(burning_seq[i] == v);
}
