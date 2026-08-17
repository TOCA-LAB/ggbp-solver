#include <optional>
#include <stdexcept>
#include <vector>

#include "ggbp/burning_simulation.hpp"
#include "ggbp/ggbp.hpp"
#include "ggbp/graph.hpp"

using std::vector;
using ggbp::simulation::BurningState;

std::optional<Vertex> best_unscheduled_source(const Graph& graph, const BurningState& state) {
    std::optional<Vertex> best_source;
    Graph::size_type best_unburned_neighbors = 0;

    for (Vertex u = 0; u < graph.order(); u++) {
        if (state.is_burned(u) || state.is_in_frontier(u))
            continue;

        const auto unburned_neighbors = graph.degree(u) - state.burned_neighbor_count(u);

        if (!best_source || unburned_neighbors >= best_unburned_neighbors) {
            best_source = u;
            best_unburned_neighbors = unburned_neighbors;
        }
    }

    return best_source;
}

Vertex choose_source(const Graph& graph, const BurningState& state) {
    if (const auto best = best_unscheduled_source(graph, state)) {
        return *best;
    }

    for (const auto u : state.frontier()) {
        if (!state.is_burned(u)) {
            return u;
        }
    }

    throw std::logic_error("unburned vertex exists but no source can be selected");
}

vector<Vertex> ggbp::heuristics::unburned_neighbors(Graph const &graph, int const r) {
    BurningState state(graph, r);
    
    vector<Vertex> burning_sequence;
    burning_sequence.reserve(graph.order());

    while (!state.complete()) {
        const auto source = choose_source(graph, state);
        burning_sequence.push_back(source);
        state.burn(source);
    }

    return burning_sequence;
}
