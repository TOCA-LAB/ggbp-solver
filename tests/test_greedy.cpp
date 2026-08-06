#include <catch2/catch_all.hpp>

#include <vector>

#include "ggbp/ggbp.hpp"
#include "ggbp/graph.hpp"
#include "support/reference_burning.hpp"

TEST_CASE("F10 tests compile against the current heuristic API",
          "[regression][build][heuristics][F10]") {
    Graph graph(1);

    const std::vector<Vertex> sequence =
        ggbp::heuristics::unburned_neighbors(graph, 1);

    REQUIRE(sequence == std::vector<Vertex>{0});
    REQUIRE(ggbp::test_support::simulate_burning(graph, 1, sequence).feasible);
}
