#include <vector>
#include <forward_list>
#include <queue>
#include <algorithm>
#include <cmath>
#include "ggbp/graph.hpp"
#include "utils/utils.hpp"

#include "ggbp/ggbp.hpp"

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
    int num_vertices = graph.order();

    vector<int> burning_sequence, to_burn;
    vector<int> burned(num_vertices, 0);
    vector<int> num_burning_neighbors(num_vertices, 0);
    vector<pair<int, double>> indexed_chromossome = indexedChromossome(chromosome);


    function<bool (pair<int, double>, pair<int, double>)> cmp = [](pair<int, double> left, pair<int, double> right) { return (left.second < right.second); };
    function<double (pair<int, double>)> val = [](pair<int, double> p) {return p.second;};
    // priority_queue<pair<int, double>, vector<pair<int, double>>, decltype(cmp)> vertex_queue(indexed_chromossome.begin(), indexed_chromossome.end(), cmp);

    ggbp::utils::bucket_sort<pair<int, double>>(indexed_chromossome, 10, cmp, val);

    auto next_to_burn = indexed_chromossome.begin();

    int n_burning = 0;
    while (n_burning < num_vertices) {
        vector<int> future_to_burn(0);
        vector<int> in_future_to_burn(num_vertices, 0);

        for (int vertex : to_burn) {
            // std::cout << vertex << std::endl;
            for (const int& neighbour : graph.neighbors(vertex)) {
                if (!burned[neighbour] && !in_future_to_burn[neighbour] && num_burning_neighbors[neighbour] >= r) {
                    future_to_burn.push_back(neighbour);
                    in_future_to_burn[neighbour] = 1;
                    n_burning++;
                }
            }
        }

        for (auto v : future_to_burn) {
            burned[v] = 1;
            for (auto u : graph.neighbors(v))
                num_burning_neighbors[u]++;
        }

        if (n_burning == num_vertices) {
            next_to_burn++;
            burning_sequence.push_back(next_to_burn->first);
        } else {
            while (next_to_burn != indexed_chromossome.end() && burned[next_to_burn->first]) {
                next_to_burn++;
            }
            if (next_to_burn != indexed_chromossome.end()) {
                burned[next_to_burn->first] = 1;
                future_to_burn.push_back(next_to_burn->first);
                burning_sequence.push_back(next_to_burn->first);
                n_burning++;

                for (auto u : graph.neighbors(next_to_burn->first))
                    num_burning_neighbors[u]++;
            }
        }
        to_burn = future_to_burn;
    }

    return { burning_sequence.size() , burning_sequence };
}

double Decoder::decode(const std::vector< double >& chromosome) const {
    return get<0>(decoder (chromosome));
}

vector<int> Decoder::getBurningSequence(const std::vector< double >& chromosome) const {
    return get<1>(decoder (chromosome));
}
}