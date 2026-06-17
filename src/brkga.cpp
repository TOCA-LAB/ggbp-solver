#include <vector>
#include <forward_list>
#include <queue>
#include <algorithm>
#include <cmath>
#include "brkga-lib/MTRand.h"
#include "brkga-lib/BRKGA.h"
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

vector<pair<int, double>> index_chromossome(const vector<double>& v)
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
    vector<pair<int, double>> indexed_chromossome = index_chromossome(chromosome);


    function<bool (pair<int, double>, pair<int, double>)> cmp = [](pair<int, double> left, pair<int, double> right) { return (left.second < right.second); };
    function<double (pair<int, double>)> val = [](pair<int, double> p) {return p.second;};

    ggbp::utils::bucket_sort<pair<int, double>>(indexed_chromossome, num_vertices / 300 + 1, cmp, val);

    auto next_to_burn = indexed_chromossome.begin();

    int n_burning = 0;
    while (n_burning < num_vertices) {
        vector<int> future_to_burn(0);
        vector<int> in_future_to_burn(num_vertices, 0);

        for (int vertex : to_burn) {
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

vector<int> Decoder::get_burning_sequence(const std::vector< double >& chromosome) const {
    return get<1>(decoder (chromosome));
}

void brkga(Graph const &graph, const int r, BRKGAParams params, vector<int> &burning_sequence, vector<int> &iters, int &first_best) {
    ggbp::brkga::Decoder decoder(graph, r);	// initialize the decoder

    // initialize the random number generator
    MTRand rng(params.seed);
    
    // initialize the BRKGA-based heuristic
    BRKGA<ggbp::brkga::Decoder, MTRand> algorithm(
        params.n,
        params.p,
        params.pe,
        params.pm,
        params.rhoe,
        decoder,
        rng,
        params.K,
        params.MAXT
    );

     unsigned generation = 0;                               // current generation
    const unsigned X_INTVL = 10;                            // exchange best individuals at every 100 generations
    const unsigned X_NUMBER = 2;                            // exchange top 2 best
    const unsigned MAX_GENS = params.MAX_GENS;              // max number of generations

    iters.reserve(MAX_GENS);

    unsigned best = graph.order();
    first_best = 0;
    
    do {
        algorithm.evolve();	// evolve the population for one generation
        int fitness = algorithm.getBestFitness();
        
        if (params.verbose == true) {
            cout << "Generation " << generation << " best solution found has objective value = "
                 << fitness << endl;
        }

        iters.push_back(fitness);

        if (fitness < best) {
            best = fitness;
            first_best = generation;
        }

        if((++generation) % X_INTVL == 0) {
            algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
        }
    } while (generation < MAX_GENS);
    
    burning_sequence = decoder.get_burning_sequence(algorithm.getBestChromosome());
}
}