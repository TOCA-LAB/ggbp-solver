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

void eligible_sources(const Graph& graph,
                      const BurningState& state,
                      vector<Vertex>& candidates) {
    candidates.clear();

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
}

struct CandidateMetrics {
    double new_frontier = 0.0;
    double weighted_deficit_reduction = 0.0;
    double residual_coverage = 0.0;
};

void project_current_frontier(
    const Graph& graph,
    const BurningState& state,
    vector<unsigned>& frontier_increments,
    vector<Vertex>& touched_vertices) {
    for (const Vertex vertex : touched_vertices) {
        frontier_increments[vertex] = 0;
    }
    touched_vertices.clear();

    const auto& burned = state.burned_vertices();
    for (const Vertex source : state.frontier()) {
        for (const Vertex neighbor : graph.neighbors(source)) {
            if (burned[neighbor] || state.is_in_frontier(neighbor)) {
                continue;
            }
            if (frontier_increments[neighbor] == 0) {
                touched_vertices.push_back(neighbor);
            }
            ++frontier_increments[neighbor];
        }
    }
}

CandidateMetrics incremental_metrics(
    const Vertex candidate,
    const Graph& graph,
    const BurningState& state,
    const vector<unsigned>& frontier_increments) {
    CandidateMetrics metrics;

    // If every remaining vertex is already in the current frontier, the
    // source is irrelevant: the whole frontier burns in this turn.
    if (state.is_in_frontier(candidate)) {
        return metrics;
    }

    const auto threshold = static_cast<unsigned>(state.threshold());
    const auto& burned = state.burned_vertices();
    const auto& burned_neighbor_counts = state.burned_neighbor_counts();

    for (const Vertex neighbor : graph.neighbors(candidate)) {
        // The candidate itself and the current frontier are marked as burned
        // before propagation, so neither can enter the next frontier.
        if (neighbor == candidate || burned[neighbor] ||
            state.is_in_frontier(neighbor)) {
            continue;
        }

        const unsigned projected_count =
            burned_neighbor_counts[neighbor] + frontier_increments[neighbor];
        if (projected_count >= threshold) {
            continue;
        }

        const unsigned deficit = threshold - projected_count;
        if (deficit == 1) {
            metrics.new_frontier += 1.0;
        }
        metrics.weighted_deficit_reduction +=
            1.0 / static_cast<double>(deficit);
        metrics.residual_coverage +=
            1.0 / static_cast<double>(projected_count + 1);
    }

    return metrics;
}

double normalize(const double value,
                 const double minimum,
                 const double maximum) {
    if (maximum == minimum) {
        return 0.0;
    }
    return (value - minimum) / (maximum - minimum);
}

void incremental_scores(
    const Graph& graph,
    const BurningState& state,
    const vector<Vertex>& candidates,
    vector<unsigned>& frontier_increments,
    vector<Vertex>& touched_vertices,
    vector<CandidateMetrics>& metrics,
    vector<double>& scores) {
    project_current_frontier(
        graph, state, frontier_increments, touched_vertices);

    metrics.clear();
    for (const Vertex candidate : candidates) {
        metrics.push_back(incremental_metrics(
            candidate, graph, state, frontier_increments));
    }

    double min_new_frontier = metrics.front().new_frontier;
    double max_new_frontier = metrics.front().new_frontier;
    double min_deficit_reduction = metrics.front().weighted_deficit_reduction;
    double max_deficit_reduction = metrics.front().weighted_deficit_reduction;
    double min_residual_coverage = metrics.front().residual_coverage;
    double max_residual_coverage = metrics.front().residual_coverage;

    for (const CandidateMetrics& candidate : metrics) {
        min_new_frontier = min(min_new_frontier, candidate.new_frontier);
        max_new_frontier = max(max_new_frontier, candidate.new_frontier);
        min_deficit_reduction =
            min(min_deficit_reduction,
                candidate.weighted_deficit_reduction);
        max_deficit_reduction =
            max(max_deficit_reduction,
                candidate.weighted_deficit_reduction);
        min_residual_coverage =
            min(min_residual_coverage, candidate.residual_coverage);
        max_residual_coverage =
            max(max_residual_coverage, candidate.residual_coverage);
    }

    constexpr double new_frontier_weight = 0.55;
    constexpr double deficit_reduction_weight = 0.30;
    constexpr double residual_coverage_weight = 0.15;

    scores.clear();
    for (const CandidateMetrics& candidate : metrics) {
        scores.push_back(
            new_frontier_weight * normalize(
                candidate.new_frontier,
                min_new_frontier,
                max_new_frontier) +
            deficit_reduction_weight * normalize(
                candidate.weighted_deficit_reduction,
                min_deficit_reduction,
                max_deficit_reduction) +
            residual_coverage_weight * normalize(
                candidate.residual_coverage,
                min_residual_coverage,
                max_residual_coverage));
    }
}

} // namespace

vector<int> construct(Graph const& graph, int r, default_random_engine& engine) {
    BurningState state(graph, r);
    vector<int> burning_sequence;
    burning_sequence.reserve(graph.order());

    // Reuse this scratch space between rounds. Only entries touched by the
    // current frontier are reset, avoiding another O(n) clear per turn.
    vector<unsigned> frontier_increments(graph.order(), 0);
    vector<Vertex> touched_vertices;
    touched_vertices.reserve(graph.order());
    vector<Vertex> candidates;
    candidates.reserve(graph.order());
    vector<Vertex> restricted_candidates;
    restricted_candidates.reserve(graph.order());
    vector<CandidateMetrics> metrics;
    metrics.reserve(graph.order());
    vector<double> scores;
    scores.reserve(graph.order());

    uniform_real_distribution<double> random_alpha(0.0, 1.0);
    const double alpha = random_alpha(engine);

    while (!state.complete()) {
        eligible_sources(graph, state, candidates);
        if (candidates.empty()) {
            throw logic_error("no eligible GRASP source");
        }

        incremental_scores(
            graph,
            state,
            candidates,
            frontier_increments,
            touched_vertices,
            metrics,
            scores);
        const auto [minimum, maximum] =
            minmax_element(scores.begin(), scores.end());
        const double cutoff = *maximum - alpha * (*maximum - *minimum);

        restricted_candidates.clear();
        for (size_t index = 0; index < candidates.size(); ++index) {
            if (scores[index] >= cutoff) {
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
