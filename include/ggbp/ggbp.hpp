#ifndef HEURISTICS_H
#define HEURISTICS_H

#include <vector>
#include <string>

#include "ggbp/graph.hpp"


namespace ggbp::heuristics {

std::vector<Vertex> greedy(const Graph &graph, const int r);
// std::vector<Graph::Vertex> unburned_neighbors(const Graph &graph, const int r);

}

namespace ggbp::ilp {

    struct {
        double ObjVal;
        double ObjBound;
        std::string msg;
    } typedef GurobiErrorInfo;


    namespace naive {
        void solve_ilp(Graph &graph, int const r, std::vector<int> &burning_sequence, GurobiErrorInfo &err_info, int verbose);
    }
    
    namespace cover {
        std::vector<int> generate_covering_tree(Graph const &g, int const r, std::vector<int> sources);

        void solve_covering_ilp(Graph const &g, int const r, int const B, std::vector<int> &burning_sequence, bool verbose);

        void ggbp_covering(Graph const &g, int const r, std::vector<int> &burning_sequence, bool verbose);    
    }
    
}

#endif
