#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "ggbp/burning_simulation.hpp"
#include "utils/utils.hpp"

#include "ggbp/ggbp.hpp"


namespace ggbp::grasp {

using namespace std;

namespace {

using ggbp::simulation::BurningState;

vector<Vertex> eligible_sources(const Graph& graph,
                                const BurningState& state) {
    vector<Vertex> candidates;
    candidates.reserve(graph.order());

    for (Vertex vertex = 0;
         static_cast<Graph::size_type>(vertex) < graph.order();
         ++vertex) {
        if (!state.is_burned(vertex) && !state.is_in_frontier(vertex)) {
            candidates.push_back(vertex);
        }
    }

    if (candidates.empty()) {
        for (const Vertex vertex : state.frontier()) {
            if (!state.is_burned(vertex)) {
                candidates.push_back(vertex);
            }
        }
    }

    return candidates;
}

Vertex greedy_completion_source(const Graph& graph,
                                const BurningState& state) {
    const auto candidates = eligible_sources(graph, state);
    if (candidates.empty()) {
        throw logic_error("no eligible GRASP source");
    }

    Vertex best = candidates.front();
    auto best_score =
        graph.degree(best) - state.burned_neighbor_count(best);

    for (const Vertex vertex : candidates) {
        const auto score =
            graph.degree(vertex) - state.burned_neighbor_count(vertex);
        if (score > best_score) {
            best = vertex;
            best_score = score;
        }
    }

    return best;
}

int completion_cost(Vertex source,
                    const Graph& graph,
                    vector<int> sequence,
                    BurningState state) {
    if (state.is_burned(source)) {
        return static_cast<int>(graph.order());
    }

    while (!state.complete()) {
        sequence.push_back(source);
        state.burn(source);

        if (!state.complete()) {
            source = greedy_completion_source(graph, state);
        }
    }

    return static_cast<int>(sequence.size());
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
    if (r < 1 || burning.size() != graph.order() ||
        num_burning_neighbors.size() != graph.order()) {
        throw invalid_argument("invalid GRASP simulation state");
    }
    if (v < 0 || static_cast<Graph::size_type>(v) >= graph.order()) {
        throw out_of_range("GRASP source is outside the graph");
    }

    BurningState state(graph, r);
    for (const Vertex source : sequence) {
        state.burn(source);
    }

    if (n_burning < 0 ||
        static_cast<Graph::size_type>(n_burning) != state.burned_count()) {
        throw invalid_argument("inconsistent GRASP burned count");
    }
    for (Vertex vertex = 0;
         static_cast<Graph::size_type>(vertex) < graph.order();
         ++vertex) {
        if (burning[vertex] != state.is_burned(vertex) ||
            num_burning_neighbors[vertex] < 0 ||
            static_cast<unsigned>(num_burning_neighbors[vertex]) !=
                state.burned_neighbor_count(vertex)) {
            throw invalid_argument("inconsistent GRASP simulation state");
        }
    }

    return completion_cost(v, graph, move(sequence), move(state));
}

vector<int> construct(Graph const& graph, int r, default_random_engine& engine) {
    BurningState state(graph, r);
    vector<int> burning_sequence;
    burning_sequence.reserve(graph.order());

    uniform_real_distribution<float> random_alpha(0, 1);
    const float alpha = random_alpha(engine);

    while (!state.complete()) {
        const auto candidates = eligible_sources(graph, state);
        if (candidates.empty()) {
            throw logic_error("no eligible GRASP source");
        }

        vector<int> costs;
        costs.reserve(candidates.size());
        for (const Vertex vertex : candidates) {
            costs.push_back(completion_cost(
                vertex, graph, burning_sequence, state));
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
        state.burn(source);
    }

    return burning_sequence;
}

namespace {

bool is_feasible(const Graph& graph,
                 const int threshold,
                 const vector<int>& sequence) {
    return ggbp::simulation::simulate_burning(
        graph, threshold, sequence).feasible;
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
