#include <vector>
#include <forward_list>
#include <queue>
#include <algorithm>
#include <cmath>
#include "ggbp/ggbp.hpp"
#include "ggbp/graph.hpp"

namespace ggbp::brkga {
using namespace std;


Decoder::Decoder(Graph const &g, int const _r) : graph(g) { 
    r = _r;
    n_vertices = g.order();
}

Decoder::~Decoder() { }

int Decoder::countBurningNeighbours(const int &v, const vector<int> &burned) const {
    int count = 0;
    for (const int& neighbour : graph.neighbors(v))
        if (burned[neighbour]) count++;
    return count;
}

vector<pair<int, double>> indexedChromossome(const vector<double>& v)
{
    vector<pair<int, double>> indexed_chromossome(v.size());
    for (unsigned i = 0; i < v.size(); ++i)
        indexed_chromossome[i] = {i, v[i]};

    return indexed_chromossome;
}

tuple<int, vector<int>> Decoder::decoder(const vector< double >& chromosome) const {
    int n_burning = 0;
    vector<int> burning_sequence, to_burn;
    vector<int> burned(n_vertices, 0);
    vector<pair<int, double>> indexed_chromossome = indexedChromossome(chromosome);

    auto cmp = [](pair<int, double> left, pair<int, double> right) { return (left.second < right.second); };
    priority_queue<pair<int, double>, vector<pair<int, double>>, decltype(cmp)> vertex_queue(indexed_chromossome.begin(), indexed_chromossome.end(), cmp);

    pair<int, double> next_to_burn = vertex_queue.top();

    while (n_burning < n_vertices) {
        vector<int> future_to_burn(0);
        vector<int> in_future_to_burn(n_vertices, 0);
        for (int vertex : to_burn) {
            // std::cout << vertex << std::endl;
            for (const int& neighbour : graph.neighbors(vertex)) {
                int burning_neighours = countBurningNeighbours(neighbour, burned);
                if (!burned[neighbour] && !in_future_to_burn[neighbour] && burning_neighours >= r) {
                    future_to_burn.push_back(neighbour);
                    in_future_to_burn[neighbour] = 1;
                    n_burning++;
                }
            }
        }

        for (auto v : future_to_burn)  
            burned[v] = 1;

        if (n_burning == n_vertices) {
            vertex_queue.pop();
            next_to_burn = vertex_queue.top();
            burning_sequence.push_back(next_to_burn.first);
        } else {
            while (vertex_queue.size() > 0 && burned[next_to_burn.first]) {
                vertex_queue.pop();
                next_to_burn = vertex_queue.top();
            }
            if (vertex_queue.size() > 0) {
                burned[next_to_burn.first] = 1;
                future_to_burn.push_back(next_to_burn.first);
                burning_sequence.push_back(next_to_burn.first);
                n_burning++;
            }
        }
        to_burn = future_to_burn;
    }

    // printVector(burning_sequence, "burning sequence");
    return { burning_sequence.size() , burning_sequence };
}

double Decoder::decode(const std::vector< double >& chromosome) const {
    return get<0>(decoder (chromosome));
}

vector<int> Decoder::getBurningSequence(const std::vector< double >& chromosome) const {
    return get<1>(decoder (chromosome));
}
}