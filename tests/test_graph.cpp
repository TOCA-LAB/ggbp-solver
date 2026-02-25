#include <catch2/catch_all.hpp>

#include "ggbp/graph.hpp"

TEST_CASE("basic construction and accessors", "[graph]") {
    Graph G(5);
    REQUIRE(G.order() == 5);
    REQUIRE(G.size() == 0);
    for (std::size_t v = 0; v < G.order(); ++v) REQUIRE(G.degree(v) == 0);
}

TEST_CASE("add and detect edges", "[graph]") {
    Graph G(4);
    G.add_edge(0,1);
    REQUIRE(G.size() == 1);
    REQUIRE(G.has_edge(0,1));
    REQUIRE(G.has_edge(1,0));
    REQUIRE(G.degree(0) == 1);
    REQUIRE(G.degree(1) == 1);

    // duplicate add doesn't increase size
    G.add_edge(0,1);
    REQUIRE(G.size() == 1);

    // different edge
    G.add_edge(2,3);
    REQUIRE(G.size() == 2);
    REQUIRE(G.has_edge(2,3));
}


TEST_CASE("graph6", "[graph]") {

    Graph G("I??@NO^Vg");

    REQUIRE(G.order() == 10);
    REQUIRE(G.size() == 17);
}



TEST_CASE("vertices iterator", "[graph]") {
    
    Graph G("I??@NO^Vg");

    int expected = 0;
    for (Vertex u : G.vertices()) {
        REQUIRE(u == expected);
        expected += 1;
    }

}
