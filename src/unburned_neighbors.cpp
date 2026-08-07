#include <optional>
#include <stdexcept>
#include <vector>

#include "ggbp/burning_simulation.hpp"
#include "ggbp/ggbp.hpp"
#include "ggbp/graph.hpp"

namespace {

std::optional<Vertex> best_unscheduled_source(
    const Graph& graph,
    const ggbp::simulation::BurningState& state) {
    std::optional<Vertex> best_source;
    Graph::size_type best_unburned_neighbors = 0;

    for (Vertex vertex = 0;
         static_cast<Graph::size_type>(vertex) < graph.order();
         ++vertex) {
        if (state.is_burned(vertex) || state.is_in_frontier(vertex)) {
            continue;
        }

        const auto unburned_neighbors =
            graph.degree(vertex) - state.burned_neighbor_count(vertex);
        if (!best_source || unburned_neighbors >= best_unburned_neighbors) {
            best_source = vertex;
            best_unburned_neighbors = unburned_neighbors;
        }
    }

    return best_source;
}

Vertex choose_source(
    const Graph& graph,
    const ggbp::simulation::BurningState& state) {
    if (const auto source = best_unscheduled_source(graph, state)) {
        return *source;
    }

    for (const Vertex vertex : state.frontier()) {
        if (!state.is_burned(vertex)) {
            return vertex;
        }
    }

    throw std::logic_error(
        "unburned vertex exists but no source can be selected");
}

} // namespace

std::vector<Vertex>
ggbp::heuristics::unburned_neighbors(const Graph& graph, const int r) {
    simulation::BurningState state(graph, r);
    std::vector<Vertex> sequence;
    sequence.reserve(graph.order());

    while (!state.complete()) {
        const Vertex source = choose_source(graph, state);
        sequence.push_back(source);
        state.burn(source);
    }

    return sequence;
}
