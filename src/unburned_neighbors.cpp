#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "ggbp/ggbp.hpp"
#include "ggbp/graph.hpp"

namespace {

struct BurningState {
    explicit BurningState(const Graph::size_type order)
        : burned(order, false),
          burned_neighbor_count(order, 0),
          in_current_frontier(order, false) {}

    std::vector<bool> burned;
    std::vector<unsigned> burned_neighbor_count;
    std::vector<Vertex> current_frontier;
    std::vector<bool> in_current_frontier;
    Graph::size_type burned_count = 0;
};

std::optional<Vertex> best_unscheduled_source(const Graph& graph,
                                               const BurningState& state) {
    std::optional<Vertex> best_source;
    Graph::size_type best_unburned_neighbors = 0;

    for (Vertex vertex = 0; static_cast<Graph::size_type>(vertex) < graph.order(); ++vertex) {
        if (state.burned[vertex] || state.in_current_frontier[vertex]) {
            continue;
        }

        const auto unburned_neighbors =
            graph.degree(vertex) - state.burned_neighbor_count[vertex];
        if (!best_source || unburned_neighbors >= best_unburned_neighbors) {
            best_source = vertex;
            best_unburned_neighbors = unburned_neighbors;
        }
    }

    return best_source;
}

Vertex choose_source(const Graph& graph, const BurningState& state) {
    if (const auto source = best_unscheduled_source(graph, state)) {
        return *source;
    }

    for (const Vertex vertex : state.current_frontier) {
        if (!state.burned[vertex]) {
            return vertex;
        }
    }

    throw std::logic_error(
        "unburned vertex exists but no source can be selected");
}

std::vector<Vertex> burning_vertices(const Vertex source,
                                     const BurningState& state) {
    std::vector<Vertex> burning_now;
    burning_now.reserve(state.current_frontier.size() + 1);
    burning_now.push_back(source);

    for (const Vertex vertex : state.current_frontier) {
        if (vertex != source && !state.burned[vertex]) {
            burning_now.push_back(vertex);
        }
    }

    return burning_now;
}

void mark_as_burned(const std::vector<Vertex>& burning_now,
                    BurningState& state) {
    for (const Vertex vertex : burning_now) {
        if (!state.burned[vertex]) {
            state.burned[vertex] = true;
            ++state.burned_count;
        }
    }
}

void build_next_frontier(const Graph& graph,
                         const int threshold,
                         const std::vector<Vertex>& burning_now,
                         BurningState& state) {
    std::vector<Vertex> next_frontier;
    std::vector<bool> in_next_frontier(graph.order(), false);

    for (const Vertex vertex : burning_now) {
        for (const Vertex neighbor : graph.neighbors(vertex)) {
            if (state.burned[neighbor]) {
                continue;
            }

            ++state.burned_neighbor_count[neighbor];
            if (state.burned_neighbor_count[neighbor] >=
                    static_cast<unsigned>(threshold) &&
                !in_next_frontier[neighbor]) {
                in_next_frontier[neighbor] = true;
                next_frontier.push_back(neighbor);
            }
        }
    }

    state.current_frontier = std::move(next_frontier);
    state.in_current_frontier = std::move(in_next_frontier);
}

void burn_round(const Graph& graph,
                const int threshold,
                const Vertex source,
                BurningState& state) {
    const auto burning_now = burning_vertices(source, state);
    mark_as_burned(burning_now, state);
    build_next_frontier(graph, threshold, burning_now, state);
}

} // namespace

std::vector<Vertex>
ggbp::heuristics::unburned_neighbors(const Graph& graph, const int r) {
    if (r < 1) {
        throw std::invalid_argument("r must be positive");
    }
    if (graph.order() == 0) {
        return {};
    }

    BurningState state(graph.order());
    std::vector<Vertex> sequence;
    sequence.reserve(graph.order());

    while (state.burned_count < graph.order()) {
        const Vertex source = choose_source(graph, state);
        sequence.push_back(source);
        burn_round(graph, r, source, state);
    }

    return sequence;
}
