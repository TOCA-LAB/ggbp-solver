#include <vector>
#include <optional>
#include "brkga-lib/MTRand.h"
#include "brkga-lib/BRKGA.h"
#include "ggbp/graph.hpp"
#include "utils/utils.hpp"
#include "ggbp/burning_simulation.hpp"

#include "ggbp/ggbp.hpp"

using ggbp::simulation::BurningState;

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

Vertex choose_source(
    const Graph& graph,
    const ggbp::simulation::BurningState& state,
    const vector<pair<int, double>> indexed_chromossome,
    vector<pair<int, double>>::iterator& candidate) {

    while (candidate != indexed_chromossome.end() &&
           state.is_burned(candidate->first)) {
        candidate++;
    }

    if (candidate != indexed_chromossome.end()) {
        const auto source = candidate->first;
        candidate++;
        return source;
    }

    for (const Vertex vertex : state.frontier()) {
        if (!state.is_burned(vertex)) {
            return vertex;
        }
    }

    throw std::logic_error("unburned vertex exists but no source can be selected");
}

tuple<int, vector<int>> Decoder::decoder(const vector< double >& chromosome) const {
    int num_vertices = graph.order();

    simulation::BurningState state(graph, r);

    vector<int> burning_sequence;
    burning_sequence.reserve(num_vertices);

    vector<pair<int, double>> indexed_chromossome = index_chromossome(chromosome);

    function<bool (pair<int, double>, pair<int, double>)> cmp = [](pair<int, double> left, pair<int, double> right) { return (left.second < right.second); };
    function<double (pair<int, double>)> val = [](pair<int, double> p) {return p.second;};

    ggbp::utils::bucket_sort<pair<int, double>>(indexed_chromossome, num_vertices / 300 + 1, cmp, val);

    auto candidate = indexed_chromossome.begin();
    while (!state.complete()) {
        const Vertex source = choose_source(graph, state, indexed_chromossome, candidate);
        burning_sequence.push_back(source);
        state.burn(source);
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
        graph.order(),
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
        
        if (params.verbose) {
            cout << "Generation " << generation + 1 << " best solution found has objective value = "
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