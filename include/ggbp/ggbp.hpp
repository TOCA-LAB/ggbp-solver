#ifndef HEURISTICS_H
#define HEURISTICS_H

#include <vector>
#include <tuple>
#include <string>
#include <random>

#include "ggbp/graph.hpp"


namespace ggbp::brkga {
    struct BRKGAParams {
        unsigned n, p, MAX_GENS;
        double pe, pm, rhoe;
        unsigned K, MAXT;
        unsigned seed;
        bool verbose;

        BRKGAParams() = default;
    };
    
    class Decoder {
    public:

        Decoder(Graph const &g, int const _r);

        ~Decoder();

        std::tuple<int, std::vector<int>> decoder(const std::vector< double >& chromosome) const;

        double decode(const std::vector< double >& chromosome) const;

        std::vector<int> get_burning_sequence(const std::vector< double >& chromosome) const;

    private:

        int r;
        Graph graph;
        int n_vertices;
    };

    void brkga(Graph const &graph, const int r, BRKGAParams params, std::vector<int> &burning_sequence, std::vector<int> &iters, int &first_best);
}

namespace ggbp::grasp {
    void grasp(
        Graph const &graph,
        int r,
        int max_it,
        std::default_random_engine &e,
        std::vector<int>& burning_sequence,
        std::vector<int> &iters,
        int &first_best,
        bool verbose
    );
}

namespace ggbp::heuristics {

std::vector<Vertex> greedy(const Graph &graph, const int r);
// std::vector<Graph::Vertex> unburned_neighbors(const Graph &graph, const int r);

}

namespace ggbp::ilp {
    struct GurobiParams {
        double mem_limit;
        double soft_mem_limit;
        double time_limit;
        bool verbose;

        GurobiParams() = default;
    };

    struct GurobiErrorInfo {
        bool error = false;
        double ObjVal;
        double ObjBound;
        std::string msg;
    };
        
    int get_upper_bound(Graph const &g, int const r, bool verbose);

    int get_lower_bound(Graph const &g, int const r, bool verbose);

    namespace naive {
        struct LPParams : ggbp::ilp::GurobiParams {
            bool linear_relax;
            bool opt;
            bool radius_cut;
            bool active_neighborhood;

            LPParams() = default;
        };
        
        void solve_ilp(Graph &g, int const r, int const B, LPParams params, std::vector<int> &burning_sequence, GurobiErrorInfo &err_info);

        void find_burning_sequence(Graph &g, int const r, int &LB, int &UB, LPParams params, std::vector<int> &burning_sequence, GurobiErrorInfo &err_info);    
    }
    
    namespace cover {    
        std::vector<int> generate_covering_tree(Graph const &g, int const r, std::vector<int> sources);

        void solve_ilp(Graph const &g, int const r, int const B, std::vector<int> &burning_sequence, ggbp::ilp::GurobiParams params, GurobiErrorInfo &err_info);

        void find_burning_sequence(Graph const &g, int const r, int LB, int UB, std::vector<int> &burning_sequence, ggbp::ilp::GurobiParams params, GurobiErrorInfo &err_info);
    }
    
    namespace prym {
        void run_bff_d(Graph &graph, std::vector<int> &burning_sequence);

        void solve_prym(Graph &graph, int B, std::vector<int> &temp_burning_sequence, std::vector<int> initial_covered_vertices);

        void run_prym(Graph &graph, std::vector<int> &burning_sequence);
    }

    namespace pap {
        void solve_ilp(Graph &graph, const int r, std::vector<int> &seed_set, ggbp::ilp::GurobiParams params);
    }
}

#endif
