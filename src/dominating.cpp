#include <forward_list>
#include <iostream>
#include <string>
#include <vector>
#include "gurobi_c++.h"
#include "ggbp/graph.hpp"
#include "ggbp/gurobi_status.hpp"

#include "ggbp/ggbp.hpp"

namespace ggbp::heuristics {
using std::vector, std::forward_list;
using std::cout, std::endl;

void min_dominating_set(Graph const &graph, vector<int> &dominating_set, ggbp::ilp::GurobiParams params, ggbp::ilp::GurobiStatus &solve_status) {
    int n_vertices = graph.order();
    
    dominating_set.clear();

    solve_status = {};
    solve_status.ObjVal = n_vertices;
    solve_status.ObjBound = 0;
    
    try {
        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, params.verbose);

        env.start();
        GRBModel model(env);

        model.set(GRB_DoubleParam_TimeLimit, params.time_limit);
        model.set(GRB_DoubleParam_SoftMemLimit, params.soft_mem_limit);
        model.set(GRB_DoubleParam_MemLimit, params.mem_limit);

        vector<GRBVar> vars_x(n_vertices);

        for (int u = 0; u < n_vertices; u++)
            vars_x[u] = model.addVar(0, 1, 0, GRB_BINARY, "x_" + std::to_string(u));
        
        for (int u = 0; u < n_vertices; u++) {
            GRBLinExpr expr = 0;
            for (auto v : graph.neighbors(u))
                expr += vars_x[v];
            model.addConstr(vars_x[u] + expr >= 1);
        }
        
        GRBLinExpr obj = 0;
        for (int u = 0; u < n_vertices; u++)
            obj += vars_x[u];
        
        model.setObjective(obj, GRB_MINIMIZE);

        model.optimize();

        auto opt_status = model.get(GRB_IntAttr_Status);
        if (opt_status == GRB_OPTIMAL) {
            dominating_set.reserve(model.get(GRB_DoubleAttr_ObjVal));
            for(int j = 0; j < n_vertices; j++) {
                if(vars_x[j].get(GRB_DoubleAttr_X) > 0.5) {
                    dominating_set.push_back(j);
                }
            }
        } else {
            report_gurobi_failure(solve_status, model, "Naive");
        }
    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;

        report_gurobi_exception(solve_status, "Naive", e.getMessage());
    } catch (...) {
        cout << "Exception during optimization" << endl;
        report_gurobi_exception(solve_status, "Naive", "unknown exception during optimization");
    }
}

vector<Vertex> dominating(Graph const &graph, int r, ggbp::ilp::GurobiParams params, ggbp::ilp::GurobiStatus &solve_status) {
    vector<int> dominating_set;

    min_dominating_set(graph, dominating_set, params, solve_status);

    if (solve_status.error) {
        return {};
    }

    int n_vertices = graph.order();
    vector<Vertex> burning_sequence;

    int n_burning = 0;
    vector<int> burning(n_vertices);
    vector<int> num_burning_neighbors(n_vertices);

    forward_list<int> to_burn;
    forward_list<int> next_to_burn;
    vector<int> in_to_burn(n_vertices, 0);

    vector<int>::iterator curr_vertex = dominating_set.begin();
    while (n_burning < n_vertices && curr_vertex != dominating_set.end()) {
        for (auto v : to_burn) {
            burning[v] = 1;
            n_burning++;
        }
        
        for (auto v : to_burn) {
            for (auto u : graph.neighbors(v)) {
                num_burning_neighbors[u]++;
                
                if (num_burning_neighbors[u] >= r && !burning[u] && !in_to_burn[u]) {
                    next_to_burn.push_front(u);
                    in_to_burn[u] = 1;
                }
            }
        }

        while (curr_vertex != dominating_set.end() && (burning[*curr_vertex] || in_to_burn[*curr_vertex])) {
            curr_vertex++;
        }
    
        if (curr_vertex == dominating_set.end())
            break;
                
        burning_sequence.push_back(*curr_vertex);
        burning[*curr_vertex] = 1;
        n_burning++;
        
        for (auto u : graph.neighbors(*curr_vertex)) {
            num_burning_neighbors[u]++;
            
            if (num_burning_neighbors[u] >= r && !burning[u] && !in_to_burn[u]) {
                next_to_burn.push_front(u);
                in_to_burn[u] = 1;
            }
        }

        curr_vertex++;

        in_to_burn = vector<int>(n_vertices, 0);
        to_burn = next_to_burn;
        next_to_burn.clear();
    }

    to_burn.clear();
    next_to_burn.clear();

    while (n_burning < n_vertices) {

        for (auto v : to_burn) {
            burning[v] = 1;
            n_burning++;
        }
        
        for (auto v : to_burn) {
            for (auto u : graph.neighbors(v)) {
                num_burning_neighbors[u]++;
                
                if (num_burning_neighbors[u] >= r && !burning[u] && !in_to_burn[u]) {
                    next_to_burn.push_front(u);
                    in_to_burn[u] = 1;
                }
            }
        }
        
        int max_vertex = 0;
        int max_neighbors = -1;

        for (int u = 0; u < n_vertices; u++) {
            if (!burning[u] && !in_to_burn[u]) {
                int new_spontaneous = 0;
                for (auto v : graph.neighbors(u)) {
                    if (num_burning_neighbors[v] == r-1 && !burning[v] && !in_to_burn[v]) {
                        new_spontaneous++;
                    }
                }

                if (new_spontaneous > max_neighbors) {
                    max_vertex = u;
                    max_neighbors = new_spontaneous;
                }
            }
        }

        burning_sequence.push_back(max_vertex);
        burning[max_vertex] = 1;
        n_burning++;

        for (auto u : graph.neighbors(max_vertex)) {
            num_burning_neighbors[u]++;
            
            if (num_burning_neighbors[u] >= r && !burning[u] && !in_to_burn[u]) {
                next_to_burn.push_front(u);
                in_to_burn[u] = 1;
            }
        }

        in_to_burn = vector<int>(n_vertices, 0);
        to_burn = next_to_burn;
        next_to_burn.clear();
    }

    return burning_sequence;
}
}
