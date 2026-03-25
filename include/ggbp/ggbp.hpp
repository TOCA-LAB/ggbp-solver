#ifndef HEURISTICS_H
#define HEURISTICS_H

#include <vector>
#include <tuple>
#include <string>

#include "ggbp/graph.hpp"


namespace ggbp::brkga {
    class Decoder {
    public:

        Decoder(Graph const &g, int const _r);

        ~Decoder();

        int countBurningNeighbours(const int &v, const std::vector<int> &burned) const;

        std::tuple<int, std::vector<int>> decoder(const std::vector< double >& chromosome) const;

        double decode(const std::vector< double >& chromosome) const;

        std::vector<int> getBurningSequence(const std::vector< double >& chromosome) const;

    private:

        int r;
        Graph graph;
        int n_vertices;
    };
}

namespace ggbp::heuristics {

std::vector<Vertex> greedy(const Graph &graph, const int r);
// std::vector<Graph::Vertex> unburned_neighbors(const Graph &graph, const int r);

}

namespace ggbp::ilp {

    struct GurobiErrorInfo {
        double ObjVal;
        double ObjBound;
        std::string msg;
    };


    namespace naive {
        void solve_ilp(Graph const &g, int const r, std::vector<int> &burning_sequence, GurobiErrorInfo &err_info, bool verbose);

        void solve_decision_ilp(Graph const &g, int const r, int const B, std::vector<int> &burning_sequence, GurobiErrorInfo &err_info, bool verbose);

        void find_burning_sequence(Graph const &g, int const r, std::vector<int> &burning_sequence, bool verbose);    
    }
    
    namespace cover {
        std::vector<int> generate_covering_tree(Graph const &g, int const r, std::vector<int> sources);

        void solve_ilp(Graph const &g, int const r, int const B, std::vector<int> &burning_sequence, bool verbose);

        void find_burning_sequence(Graph const &g, int const r, std::vector<int> &burning_sequence, bool verbose);    
    }
    
    namespace prym {
        void run_bff_d(Graph &graph, std::vector<int> &burning_sequence);

        void solve_prym(Graph &graph, int B, std::vector<int> &temp_burning_sequence, std::vector<int> initial_covered_vertices);

        void run_prym(Graph &graph, std::vector<int> &burning_sequence);
    }
}

#endif
