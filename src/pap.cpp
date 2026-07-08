#include "gurobi_c++.h"
#include <vector>
#include <iostream>
#include "ggbp/graph.hpp"

#include "ggbp/ggbp.hpp"

namespace ggbp::ilp::pap {
using std::vector;
using std::cout, std::endl;

void solve_ilp(Graph &graph, const int r, std::vector<int> &seed_set, ggbp::ilp::GurobiParams params, ggbp::ilp::GurobiErrorInfo &err_info) {
    int const n_vertices = graph.order();

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

        int opt_status = model.get(GRB_IntAttr_Status);

        seed_set.clear();

        if (opt_status == GRB_OPTIMAL) {
            for(int u = 0; u < n_vertices; u++) {
                if (s[u][0].get(GRB_DoubleAttr_X) > 0.5) {
                    seed_set.push_back(u);
                }
            }
        } else if (opt_status == GRB_TIME_LIMIT) {
            err_info.error = true;
            
            if (model.get(GRB_DoubleAttr_ObjVal) < n_vertices)
                err_info.ObjVal = model.get(GRB_DoubleAttr_ObjVal);
            else
                err_info.ObjVal = n_vertices;
            
            err_info.ObjBound = model.get(GRB_DoubleAttr_ObjBound);
            err_info.msg = "time limit";
        } else if (opt_status == GRB_MEM_LIMIT) {
            err_info.error = true;

            if (model.get(GRB_DoubleAttr_ObjVal) < n_vertices)
                err_info.ObjVal = model.get(GRB_DoubleAttr_ObjVal);
            else
                err_info.ObjVal = n_vertices;

            err_info.ObjBound = model.get(GRB_DoubleAttr_ObjBound);
            err_info.msg = "memory limit";
        }

    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;

        if (e.getErrorCode() == GRB_ERROR_OUT_OF_MEMORY) {
            err_info.error = true;
            err_info.ObjVal = n_vertices;
            err_info.ObjBound = 0;
            err_info.msg = "out of memory";
        }
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }
}

}