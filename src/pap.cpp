#include "gurobi_c++.h"
#include <vector>
#include <iostream>
#include "ggbp/graph.hpp"
#include "ggbp/gurobi_status.hpp"

#include "ggbp/ggbp.hpp"

namespace ggbp::ilp::pap {
using std::vector;
using std::cout, std::endl;

namespace {

bool is_valid_uniform_pap_seed_set(
    const Graph& graph,
    const int threshold,
    const vector<int>& seed_set) {
    if (threshold <= 0) {
        return false;
    }

    vector<bool> seeds(graph.order(), false);
    vector<bool> spreaders(graph.order(), false);
    for (const int seed : seed_set) {
        if (seed < 0 ||
            static_cast<Graph::size_type>(seed) >= graph.order()) {
            return false;
        }
        seeds[seed] = true;
        spreaders[seed] = true;
    }

    for (Graph::size_type turn = 0; turn < graph.order(); ++turn) {
        vector<bool> next = spreaders;
        bool changed = false;
        for (Vertex vertex = 0;
             static_cast<Graph::size_type>(vertex) < graph.order();
             ++vertex) {
            if (spreaders[vertex]) {
                continue;
            }

            int active_neighbors = 0;
            for (const Vertex neighbor : graph.neighbors(vertex)) {
                active_neighbors += spreaders[neighbor] ? 1 : 0;
            }
            if (active_neighbors >= threshold) {
                next[vertex] = true;
                changed = true;
            }
        }
        spreaders = std::move(next);
        if (!changed) {
            break;
        }
    }

    for (Vertex vertex = 0;
         static_cast<Graph::size_type>(vertex) < graph.order();
         ++vertex) {
        if (seeds[vertex]) {
            continue;
        }
        bool aware = false;
        for (const Vertex neighbor : graph.neighbors(vertex)) {
            aware = aware || spreaders[neighbor];
        }
        if (!aware) {
            return false;
        }
    }
    return true;
}

} // namespace

void solve_ilp(Graph &graph, const int r, std::vector<int> &seed_set, ggbp::ilp::GurobiParams params, ggbp::ilp::GurobiErrorInfo &err_info) {
    int const n_vertices = graph.order();
    seed_set.clear();
    err_info = {};

    try {        
        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, params.verbose);
        env.start();
        GRBModel model(env);

        model.set(GRB_DoubleParam_TimeLimit, params.time_limit);
        model.set(GRB_DoubleParam_SoftMemLimit, params.soft_mem_limit);
        model.set(GRB_DoubleParam_MemLimit, params.mem_limit);

        // --- Variáveis ---
        // s[v, t]: 1 se o vértice v é disseminador no turno t
        vector<vector<GRBVar>> s(n_vertices, vector<GRBVar>(n_vertices + 1));

        for (int u = 0; u < n_vertices; u++) {
            for (int t = 0; t < n_vertices + 1; t++) {
                s[u][t] = model.addVar(0, 1, 0, GRB_BINARY, "s_" + std::to_string(u) + "_" + std::to_string(t));
            }
        }

        // --- Restrições ---
        // 1. Para u ser disseminador no turno t, ou está no conjunto semente,
        // ou tem r vizinhos disseminadores no turno t-1
        for (int u = 0; u < n_vertices; u++) {
            for (int t = 1; t < n_vertices + 1; t++) {
                GRBLinExpr expr = 0;
                
                for (auto v : graph.neighbors(u)) {
                    expr += s[v][t-1];
                }

                model.addConstr(expr - r * (s[u][t] - s[u][0]) >= 0, "Spreader_" + std::to_string(u) + "_" + std::to_string(t));
            }
        }

        // 2. Ao final, todo vértice deve ser conhecedor
        for (int u = 0; u < n_vertices; u++) {
            GRBLinExpr expr = 0;
            for (auto v : graph.neighbors(u)) {
                expr += s[v][n_vertices - 1];
            }

            model.addConstr(s[u][0] + expr >= 1, "Aware_" + std::to_string(u));
        }

        // 3. Restrição adicional
        for (int u = 0; u < n_vertices; u++) {
            for (int t = 1; t < n_vertices + 1; t++) {
                GRBLinExpr expr = 0;

                for (auto v : graph.neighbors(u)) {
                    expr += s[v][t];
                }
                
                for (int v = 0; v < n_vertices; v++) {
                    expr += s[v][t] - s[v][t-1];
                }
                model.addConstr(s[u][t] + expr >= 1, "Adicional_" + std::to_string(u) + "_" + std::to_string(t));
            }
        }

        // Minimiza disseminadores iniciais (conjunto semente)
        GRBLinExpr obj = 0;
        for (int u = 0; u < n_vertices; u++) {
            obj += s[u][0];
        }
        
        model.setObjective(obj, GRB_MINIMIZE);

        model.optimize();

        const auto solve_status = inspect_gurobi_model(model);
        record_gurobi_status(err_info, solve_status);

        if (solve_status.status == GurobiSolveStatus::Optimal ||
            solve_status.status == GurobiSolveStatus::Feasible) {
            for(int u = 0; u < n_vertices; u++) {
                if (s[u][0].get(GRB_DoubleAttr_X) > 0.5) {
                    seed_set.push_back(u);
                }
            }

            if (!is_valid_uniform_pap_seed_set(graph, r, seed_set)) {
                seed_set.clear();
                err_info.error = true;
                err_info.status = GurobiSolveStatus::Error;
                err_info.msg = "PAP returned an invalid seed-set incumbent";
                return;
            }

            if (solve_status.status == GurobiSolveStatus::Feasible) {
                err_info.ObjVal = model.get(GRB_DoubleAttr_ObjVal);
                err_info.ObjBound = model.get(GRB_DoubleAttr_ObjBound);
                err_info.error = true;
                err_info.msg =
                    "PAP retained a non-optimal incumbent at " +
                    std::string(gurobi_status_name(
                        solve_status.native_status));
            }
        } else {
            report_gurobi_failure(err_info, solve_status, "PAP");
        }

    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;

        report_gurobi_exception(err_info, "PAP", e.getMessage());
        err_info.ObjVal = n_vertices;
        err_info.ObjBound = 0;
    } catch (...) {
        cout << "Exception during optimization" << endl;
        report_gurobi_exception(
            err_info, "PAP", "unknown exception during optimization");
        err_info.ObjVal = graph.order();
        err_info.ObjBound = 0;
    }
}

}
