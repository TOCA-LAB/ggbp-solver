#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <random>
#include <forward_list>

#include "utils/utils.hpp"

#include "ggbp/ggbp.hpp"


namespace ggbp::grasp {

using namespace std;

// Use the unburned_neighbors heuristic to estimate how good a vertex is
int calculate_cost(int v, Graph const &graph, int r, vector<int> sequence, vector<bool> burning, int n_burning, vector<int> num_burning_neighbors) {
    int n_vertices = graph.order();
    
    if (burning[v])
        return n_vertices;

    forward_list<int> to_burn;

    bool first = true;
    int max_vertex = v;
    

    while (n_burning < n_vertices) {
        forward_list<int> next_to_burn;
        vector<int> in_to_burn(n_vertices, 0);

        int max_neighbors = 0;

        if (first) {
            first = false;
        } else {
            for (int u = 0; u < n_vertices; u++) {
                if (graph.degree(u) - num_burning_neighbors[u] > max_neighbors && !burning[u] && !in_to_burn[u]) {
                    max_vertex = u;
                    max_neighbors = graph.degree(u) - num_burning_neighbors[u];
                }
            }
        }

        sequence.push_back(max_vertex);
        to_burn.push_front(max_vertex);
        in_to_burn[max_vertex] = 1;

        if (n_burning == n_vertices)
            break;
    
        for (auto v : to_burn) {
            burning[v] = 1;
            n_burning++;
        }
        
        for (auto v : to_burn) {
            for (auto u : graph.neighbors(v)) {
                num_burning_neighbors[u]++;
                
                if (num_burning_neighbors[u] >= r && !burning[u] && !in_to_burn[u]) {
                    next_to_burn.push_front(u);
                    in_to_burn[u] = 1;
                }
            }
        }

        if (n_burning == n_vertices)
            sequence.push_back(n_vertices - 1);
        
        to_burn = next_to_burn;
        next_to_burn.clear();
    }

    return sequence.size();
}

vector<int> construct(Graph const &graph, int r, default_random_engine &e) {
    int n_vertices = graph.order(), n_burning = 0;

    vector<int> burning_sequence, cost(n_vertices), num_burning_neighbors(n_vertices);
    forward_list<int> to_burn, next_to_burn;
    vector<bool> burning(n_vertices, 0), in_to_burn(n_vertices, 0);

    uniform_real_distribution<float> f(0,1);
    float a = f(e);
    
    while (n_burning < n_vertices) {
        cost[0] = calculate_cost(0, graph, r, burning_sequence, burning, n_burning, num_burning_neighbors); 

        int c_min = cost[0], c_max = cost[0];
        for (int i = 1; i < n_vertices; i++) {
            cost[i] = calculate_cost(i, graph, r, burning_sequence, burning, n_burning, num_burning_neighbors);

            if (cost[i] < c_min)
                c_min = cost[i];
            if (cost[i] > c_max)
                c_max = cost[i];
        }

        float threshold = c_min + a * (c_max - c_min);

        // Populate the restricted candidate list using a quality-based scheme
        vector<int> RCL;
        RCL.reserve(n_vertices/10);

        for (int i = 0; i < n_vertices; i++)
            if (cost[i] <= threshold)
                RCL.push_back(i);

        RCL.shrink_to_fit();

        // Choose random unburned vertex from RCL
        uniform_int_distribution<> u(0, RCL.size()-1);
        int current;
        do {
            current = RCL[u(e)];
        } while (burning[current]);
    
        // Burn the choosen vertex and propagate the fire
        burning_sequence.push_back(current);
        if (!in_to_burn[current]) {
            to_burn.push_front(current);
            in_to_burn[current] = true;
        }

        for (auto v : to_burn) {
            burning[v] = 1;
            n_burning++;
        }

        for (auto v : to_burn) {
            for (auto u : graph.neighbors(v)) {
                num_burning_neighbors[u]++;
                
                if (num_burning_neighbors[u] >= r && !burning[u] && !in_to_burn[u]) {
                    next_to_burn.push_front(u);
                    in_to_burn[u] = true;
                }
            }
        }

        to_burn = next_to_burn;
        next_to_burn.clear();
    }
    
    return burning_sequence;
}

void local_search(vector<int>& burning_sequence, Graph const &graph, int r) {
    int n_vertices = graph.order();

    bool changed = true;
    
    while(changed) {
        changed = false;

        vector<int> temp_burning_sequence;

        vector<int> best_burning_sequence = burning_sequence;
        vector<int> best_history;

        ggbp::utils::burning_history(graph, r, best_burning_sequence, best_history);
        
        for (int i = 0; i < burning_sequence.size(); i++) {
            temp_burning_sequence = burning_sequence;
            temp_burning_sequence.erase(temp_burning_sequence.begin() + i);

            vector<int> history;
            bool valid_sequence = ggbp::utils::burning_history(graph, r, temp_burning_sequence, history);

            if (valid_sequence) {
                if (temp_burning_sequence.size() < best_burning_sequence.size() && history > best_history) {
                    swap(best_burning_sequence, temp_burning_sequence);
                    swap(best_history, history);
                }
            }
        }

        if (best_burning_sequence.size() < burning_sequence.size()) {
            swap(burning_sequence, best_burning_sequence);
            changed = true;
        }
    }

    vector<bool> in_sequence(n_vertices, false);

    for (auto u : burning_sequence) {
        in_sequence[u] = true;
    }

    int k = (n_vertices - burning_sequence.size()) / 10;

    auto cmp = [&graph](int left, int right) { return (graph.degree(left) > graph.degree(right)); };

    vector<int> vertex_queue;
    vertex_queue.reserve(k);

    for (int u = 0; u < n_vertices; u++) {
        if (!in_sequence[u]) {
            if (vertex_queue.size() < k) {
                vertex_queue.push_back(u);
                push_heap(vertex_queue.begin(), vertex_queue.end());
            } else {
                pop_heap(vertex_queue.begin(), vertex_queue.end());
                vertex_queue.pop_back();

                vertex_queue.push_back(u);
                push_heap(vertex_queue.begin(), vertex_queue.end());
            }
        }
    }

    for (auto u : vertex_queue) {
        vector<int> temp_burning_sequence;

        vector<int> best_burning_sequence = burning_sequence;
        int best_turns = burning_sequence.size();
        
        for (int i = 0; i < burning_sequence.size(); i++) {
            temp_burning_sequence = burning_sequence;
            temp_burning_sequence[i] = u;

            int turns = ggbp::utils::burning_turns(graph, r, temp_burning_sequence);

            if (turns != 0 && turns < best_turns) {
                swap(best_burning_sequence, temp_burning_sequence);
                best_turns = turns;
            }
        }

        best_burning_sequence.resize(best_turns);

        if (best_burning_sequence.size() < burning_sequence.size()) {
            swap(burning_sequence, best_burning_sequence);
        }
    }

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
