#include <vector>
#include <forward_list>
#include <algorithm>
#include <cmath>
#include <functional>

#include "ggbp/graph.hpp"

#include "utils/utils.hpp"

namespace ggbp::utils {

using std::vector, std::forward_list;

bool is_burning_sequence(Graph const &g, int r, vector<int> &burning_sequence, Vertex &unburned) {
    int num_vertices = g.order();

    vector<bool> in_sequence(num_vertices, false);
    vector<bool> burned(num_vertices, false);
    vector<int> num_burned_neighbors(num_vertices, 0);


    forward_list<int> to_burn;
    vector<bool> in_to_burn(num_vertices, false);
    forward_list<int> next_to_burn;

    for (auto u : burning_sequence) {
        if (!burned[u]) {
            burned[u] = true;

            for (auto w : g.neighbors(u)) {
                num_burned_neighbors[w]++;

                if (!in_to_burn[w] && !burned[w] && num_burned_neighbors[w] >= r) {
                    next_to_burn.push_front(w);
                    in_to_burn[w] = true;
                }
            }
        }

        for (auto v : to_burn) {
            if (!burned[v]) {
                burned[v] = true;
                for (auto w : g.neighbors(v)) {
                    num_burned_neighbors[w]++;

                    if (!in_to_burn[w] && !burned[w] && num_burned_neighbors[w] >= r) {
                        next_to_burn.push_front(w);
                        in_to_burn[w] = true;
                    }
                }
            }
        }
        
        to_burn = next_to_burn;
        next_to_burn.clear();
        in_to_burn.assign(num_vertices, false);
    }

    for (int u = 0; u < num_vertices; u++) {
        if (!burned[u]) {
            unburned = u;
            return false;
        }
    }

    return true;
}

int burning_turns(Graph const &g, int r, vector<int> const &burning_sequence) {
    int num_vertices = g.order();

    int num_burned = 0;
    vector<bool> burned(num_vertices, false);
    vector<int> num_burned_neighbors(num_vertices, 0);


    forward_list<int> to_burn;
    vector<bool> in_to_burn(num_vertices, false);
    forward_list<int> next_to_burn;

    for (int i = 0; i < burning_sequence.size(); i++) {
        auto u = burning_sequence[i];

        if (!burned[u]) {
            burned[u] = true;
            num_burned++;

            for (auto w : g.neighbors(u)) {
                num_burned_neighbors[w]++;

                if (!in_to_burn[w] && !burned[w] && num_burned_neighbors[w] >= r) {
                    next_to_burn.push_front(w);
                    in_to_burn[w] = true;
                }
            }
        }

        for (auto v : to_burn) {
            if (!burned[v]) {
                burned[v] = true;
                num_burned++;

                for (auto w : g.neighbors(v)) {
                    num_burned_neighbors[w]++;

                    if (!in_to_burn[w] && !burned[w] && num_burned_neighbors[w] >= r) {
                        next_to_burn.push_front(w);
                        in_to_burn[w] = true;
                    }
                }
            }
        }
    
        if (num_burned >= num_vertices) {
            return i + 1;
        }

        to_burn = next_to_burn;
        next_to_burn.clear();
        in_to_burn.assign(num_vertices, false);
    }

    for (int u = 0; u < num_vertices; u++) {
        if (!burned[u]) {
            return 0;
        }
    }

    return burning_sequence.size();
}

bool burning_history(Graph const &g, int r, vector<int> const &burning_sequence, vector<int> &history) {
    int num_vertices = g.order();
    int B = burning_sequence.size();

    vector<bool> in_sequence(num_vertices, false);
    vector<bool> burned(num_vertices, false);
    vector<int> num_burned_neighbors(num_vertices, 0);

    history.assign(B,0);

    forward_list<int> to_burn;
    vector<bool> in_to_burn(num_vertices, false);
    forward_list<int> next_to_burn;

    for (int i = 0; i < B; i++) {
        auto u = burning_sequence[i];
        
        if (!burned[u]) {
            burned[u] = true;
            history[i]++;


            for (auto w : g.neighbors(u)) {
                num_burned_neighbors[w]++;

                if (!in_to_burn[w] && !burned[w] && num_burned_neighbors[w] >= r) {
                    next_to_burn.push_front(w);
                    in_to_burn[w] = true;
                }
            }
        }

        for (auto v : to_burn) {
            if (!burned[v]) {
                burned[v] = true;
                history[i]++;

                for (auto w : g.neighbors(v)) {
                    num_burned_neighbors[w]++;

                    if (!in_to_burn[w] && !burned[w] && num_burned_neighbors[w] >= r) {
                        next_to_burn.push_front(w);
                        in_to_burn[w] = true;
                    }
                }
            }
        }
        
        to_burn = next_to_burn;
        next_to_burn.clear();
        in_to_burn.assign(num_vertices, false);
    }

    for (int u = 0; u < num_vertices; u++) {
        if (!burned[u])
            return false;
    }

    return true;
}

}