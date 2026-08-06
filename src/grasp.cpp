#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "utils/utils.hpp"

#include "ggbp/ggbp.hpp"


namespace ggbp::grasp {

using namespace std;

namespace {

struct ConstructionState {
    explicit ConstructionState(const Graph::size_type order)
        : burned(order, false),
          burned_neighbor_count(order, 0),
          in_current_frontier(order, false) {}

    vector<bool> burned;
    vector<int> burned_neighbor_count;
    vector<Vertex> current_frontier;
    vector<bool> in_current_frontier;
    int burned_count = 0;
};

void rebuild_frontier_from_counts(const Graph& graph,
                                  const int threshold,
                                  ConstructionState& state) {
    state.current_frontier.clear();
    state.in_current_frontier.assign(graph.order(), false);

    for (Vertex vertex = 0;
         static_cast<Graph::size_type>(vertex) < graph.order();
         ++vertex) {
        if (!state.burned[vertex] &&
            state.burned_neighbor_count[vertex] >= threshold) {
            state.current_frontier.push_back(vertex);
            state.in_current_frontier[vertex] = true;
        }
    }
}

vector<Vertex> eligible_sources(const Graph& graph,
                                const ConstructionState& state) {
    vector<Vertex> candidates;
    candidates.reserve(graph.order());

    for (Vertex vertex = 0;
         static_cast<Graph::size_type>(vertex) < graph.order();
         ++vertex) {
        if (!state.burned[vertex] && !state.in_current_frontier[vertex]) {
            candidates.push_back(vertex);
        }
    }

    if (candidates.empty()) {
        for (const Vertex vertex : state.current_frontier) {
            if (!state.burned[vertex]) {
                candidates.push_back(vertex);
            }
        }
    }

    return candidates;
}

vector<Vertex> burning_vertices(const Vertex source,
                                const ConstructionState& state) {
    vector<Vertex> burning_now;
    burning_now.reserve(state.current_frontier.size() + 1);
    burning_now.push_back(source);

    for (const Vertex vertex : state.current_frontier) {
        if (vertex != source && !state.burned[vertex]) {
            burning_now.push_back(vertex);
        }
    }

    return burning_now;
}

void burn_round(const Graph& graph,
                const int threshold,
                const Vertex source,
                ConstructionState& state) {
    const auto burning_now = burning_vertices(source, state);

    for (const Vertex vertex : burning_now) {
        if (!state.burned[vertex]) {
            state.burned[vertex] = true;
            ++state.burned_count;
        }
    }

    vector<Vertex> next_frontier;
    vector<bool> in_next_frontier(graph.order(), false);

    for (const Vertex vertex : burning_now) {
        for (const Vertex neighbor : graph.neighbors(vertex)) {
            if (state.burned[neighbor]) {
                continue;
            }

            ++state.burned_neighbor_count[neighbor];
            if (state.burned_neighbor_count[neighbor] >= threshold &&
                !in_next_frontier[neighbor]) {
                next_frontier.push_back(neighbor);
                in_next_frontier[neighbor] = true;
            }
        }
    }

    state.current_frontier = move(next_frontier);
    state.in_current_frontier = move(in_next_frontier);
}

Vertex greedy_completion_source(const Graph& graph,
                                const ConstructionState& state) {
    const auto candidates = eligible_sources(graph, state);
    if (candidates.empty()) {
        throw logic_error("no eligible GRASP source");
    }

    Vertex best = candidates.front();
    auto best_score =
        graph.degree(best) - state.burned_neighbor_count[best];

    for (const Vertex vertex : candidates) {
        const auto score =
            graph.degree(vertex) - state.burned_neighbor_count[vertex];
        if (score > best_score) {
            best = vertex;
            best_score = score;
        }
    }

    return best;
}

} // namespace

// Use the Unburned Neighbors rule to estimate the completed sequence length.
int calculate_cost(int v,
                   Graph const& graph,
                   int r,
                   vector<int> sequence,
                   vector<bool> burning,
                   int n_burning,
                   vector<int> num_burning_neighbors) {
    const int n_vertices = graph.order();
    if (r < 1 || burning.size() != graph.order() ||
        num_burning_neighbors.size() != graph.order()) {
        throw invalid_argument("invalid GRASP simulation state");
    }
    if (v < 0 || v >= n_vertices) {
        throw out_of_range("GRASP source is outside the graph");
    }
    if (burning[v]) {
        return n_vertices;
    }

    ConstructionState state(graph.order());
    state.burned = move(burning);
    state.burned_neighbor_count = move(num_burning_neighbors);
    state.burned_count = n_burning;
    rebuild_frontier_from_counts(graph, r, state);

    Vertex source = v;
    while (state.burned_count < n_vertices) {
        sequence.push_back(source);
        burn_round(graph, r, source, state);

        if (state.burned_count < n_vertices) {
            source = greedy_completion_source(graph, state);
        }
    }

    return sequence.size();
}

vector<int> construct(Graph const& graph, int r, default_random_engine& engine) {
    ConstructionState state(graph.order());
    vector<int> burning_sequence;
    burning_sequence.reserve(graph.order());

    uniform_real_distribution<float> random_alpha(0, 1);
    const float alpha = random_alpha(engine);

    while (state.burned_count < static_cast<int>(graph.order())) {
        const auto candidates = eligible_sources(graph, state);
        if (candidates.empty()) {
            throw logic_error("no eligible GRASP source");
        }

        vector<int> costs;
        costs.reserve(candidates.size());
        for (const Vertex vertex : candidates) {
            costs.push_back(calculate_cost(vertex,
                                           graph,
                                           r,
                                           burning_sequence,
                                           state.burned,
                                           state.burned_count,
                                           state.burned_neighbor_count));
        }

        const auto [minimum, maximum] = minmax_element(costs.begin(), costs.end());
        const float cutoff = *minimum + alpha * (*maximum - *minimum);

        vector<Vertex> restricted_candidates;
        restricted_candidates.reserve(candidates.size());
        for (size_t index = 0; index < candidates.size(); ++index) {
            if (costs[index] <= cutoff) {
                restricted_candidates.push_back(candidates[index]);
            }
        }

        uniform_int_distribution<size_t> choose(
            0, restricted_candidates.size() - 1);
        const Vertex source = restricted_candidates[choose(engine)];

        burning_sequence.push_back(source);
        burn_round(graph, r, source, state);
    }

    return burning_sequence;
}

namespace {

bool is_feasible(const Graph& graph,
                 const int threshold,
                 const vector<int>& sequence) {
    vector<int> history;
    return ggbp::utils::burning_history(graph, threshold, sequence, history);
}

void remove_redundant_sources(const Graph& graph,
                              const int threshold,
                              vector<int>& sequence) {
    bool removed_source;
    do {
        removed_source = false;

        for (size_t index = 0; index < sequence.size(); ++index) {
            auto candidate = sequence;
            candidate.erase(candidate.begin() + index);

            if (is_feasible(graph, threshold, candidate)) {
                sequence = move(candidate);
                removed_source = true;
                break;
            }
        }
    } while (removed_source);
}

size_t replacement_candidate_limit(const Graph& graph,
                                   const vector<int>& sequence) {
    if (sequence.size() >= graph.order()) {
        return 0;
    }

    const size_t vertices_outside_sequence = graph.order() - sequence.size();
    const size_t requested = vertices_outside_sequence / 10 +
                             (vertices_outside_sequence < 10 ? 2 : 0);
    return min(requested, vertices_outside_sequence);
}

vector<Vertex> select_replacement_candidates(const Graph& graph,
                                             const vector<int>& sequence) {
    const size_t limit = replacement_candidate_limit(graph, sequence);
    if (limit == 0) {
        return {};
    }

    vector<bool> in_sequence(graph.order(), false);
    for (const Vertex source : sequence) {
        in_sequence[source] = true;
    }

    auto cmp = [&graph](const Vertex left, const Vertex right) {
        return graph.degree(left) > graph.degree(right);
    };

    vector<Vertex> candidates;
    candidates.reserve(limit);

    for (Vertex vertex = 0;
         static_cast<Graph::size_type>(vertex) < graph.order();
         ++vertex) {
        if (in_sequence[vertex]) {
            continue;
        }

        if (candidates.size() < limit) {
            candidates.push_back(vertex);
            push_heap(candidates.begin(), candidates.end(), cmp);
            continue;
        }

        if (graph.degree(vertex) > graph.degree(candidates.front())) {
            pop_heap(candidates.begin(), candidates.end(), cmp);
            candidates.back() = vertex;
            push_heap(candidates.begin(), candidates.end(), cmp);
        }
    }

    return candidates;
}

bool find_shortening_replacement(const Graph& graph,
                                 const int threshold,
                                 const vector<int>& sequence,
                                 const Vertex replacement,
                                 vector<int>& improved_sequence) {
    size_t best_turns = sequence.size();
    bool found_improvement = false;

    for (size_t index = 0; index < sequence.size(); ++index) {
        auto candidate = sequence;
        candidate[index] = replacement;

        const int turns =
            ggbp::utils::burning_turns(graph, threshold, candidate);
        if (turns > 0 && static_cast<size_t>(turns) < best_turns) {
            best_turns = turns;
            candidate.resize(best_turns);
            improved_sequence = move(candidate);
            found_improvement = true;
        }
    }

    return found_improvement;
}

void apply_replacement_improvements(const Graph& graph,
                                    const int threshold,
                                    vector<int>& sequence) {
    const auto replacements = select_replacement_candidates(graph, sequence);

    for (const Vertex replacement : replacements) {
        vector<int> improved_sequence;
        if (find_shortening_replacement(
                graph, threshold, sequence, replacement, improved_sequence)) {
            sequence = move(improved_sequence);
        }
    }
}

} // namespace

void local_search(vector<int>& burning_sequence, Graph const& graph, int r) {
    remove_redundant_sources(graph, r, burning_sequence);
    apply_replacement_improvements(graph, r, burning_sequence);
}


void grasp(
    Graph const &graph,
    int r,
    int max_it,
    default_random_engine &e,
    vector<int> &burning_sequence,
    vector<int> &iters,
    int &first_best,
    bool verbose)
{
    vector<int> temp_burning_sequence;
    iters.reserve(max_it);
    first_best = 0;
    int best = graph.order();

    for (int i = 0; i < max_it; i++) {
        if (verbose)
            cout << "Iteration " << i+1 << endl;
        
        temp_burning_sequence = construct(graph, r, e);

        local_search(temp_burning_sequence, graph, r);

        if (burning_sequence.size() == 0)
            swap(temp_burning_sequence, burning_sequence);
        else if (temp_burning_sequence.size() < burning_sequence.size())
            swap(temp_burning_sequence, burning_sequence);
        
        if (verbose)
            cout << "Best: " << burning_sequence.size() << endl;

        iters.push_back(burning_sequence.size());

        if (iters[i] < best) {
            best = iters[i];
            first_best = i;
        }
    }
}

}
