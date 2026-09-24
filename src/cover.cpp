#include <vector>
#include <queue>
#include <forward_list>
#include <iostream>

#include <gurobi_c++.h>
#include "ggbp/graph.hpp"
#include "utils/combinations.hpp"
#include "utils/contract.hpp"
#include "utils/utils.hpp"

#include "ggbp/ggbp.hpp"

namespace ggbp::ilp::cover {
using namespace std;

vector<int> generate_covering_tree(Graph const &g, int const r, vector<Vertex> sources) {
    vector<int> depth(g.order(), -1);
    vector<Vertex> tree_neighbors(g.order(), 0);

    queue<int> q;
    for (auto src : sources) {
        depth[src] = 0;
        q.push(src);
    }

    while (!q.empty()) {
        Vertex u = q.front();
        q.pop();

        auto neighbors = g.neighbors(u);
        for (Vertex v : neighbors) {
            if (depth[v] != -1)
                continue;
        
            tree_neighbors[v]++;

            if (tree_neighbors[v] >= r) {
                depth[v] = depth[u] + 1;
                q.push(v);
            }
        }
    }

    return depth;
}

void solve_ilp(Graph const &g, int const r, int const B, ggbp::ilp::GurobiParams params, vector<int> &burning_sequence, GurobiStatus &solve_status) {
    int num_vertices = g.order();
    vector<int> vertices = g.vertices();

    burning_sequence.clear();

    solve_status = {};
    solve_status.ObjVal = B;
    solve_status.ObjBound = 0;

    try {
        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, params.verbose);
        env.start();

        GRBModel model(env);
        
        model.set(GRB_DoubleParam_TimeLimit, params.time_limit);
        model.set(GRB_DoubleParam_SoftMemLimit, params.soft_mem_limit);
        model.set(GRB_DoubleParam_MemLimit, params.mem_limit);

        vector<vector<GRBVar>> vars_x(num_vertices, vector<GRBVar>(B));
        vector<vector<GRBVar>> vars_y(num_vertices, vector<GRBVar>(B));
        vector<vector<GRBVar>> vars_f;

        for (Vertex u = 0; u < num_vertices; u++) {
            for (int i = 0; i < B; i++) {
                vars_x[u][i] = model.addVar(0, 1, 0, GRB_BINARY, "x_" + to_string(u) + "_" + to_string(i));
            }
        }

        for (Vertex u = 0; u < num_vertices; u++) {
            for (int i = 0; i < B; i++) {
                vars_y[u][i] = model.addVar(0, 1, 0, GRB_BINARY, "y_" + to_string(u) + "_" + to_string(i));
            }
        }

        for (int i = 0; i < B; i++) {
            GRBLinExpr expr = 0;
            for (int u = 0; u < num_vertices; u++) {
                expr += vars_x[u][i];
            }
            model.addConstr(expr == 1);
        }

        for (int u = 0; u < num_vertices; u++) {
            GRBLinExpr expr = 0;
            for (int i = 0; i < B; i++) {
                expr += vars_x[u][i];
            }
            model.addConstr(expr <= 1);
        }

        for (int u = 0; u < num_vertices; u++) {
            model.addConstr(vars_y[u][B-1] == 1);
        }

        auto combs = ggbp::utils::Combinations(vertices, r);

        unsigned long long s = 0;
        for (auto comb : combs) {
            vars_f.push_back(vector<GRBVar>(B));
            for (int i = 0; i < B; i++) {
                vars_f[s][i] = model.addVar(0, 1, 0, GRB_BINARY, "F_" + to_string(s) + "_" + to_string(i));
            }

            for (int i = 1; i < B; i++) {
                GRBLinExpr expr = 0;
                for (auto u : comb) {
                    expr += vars_y[u][i-1];
                }
                model.addConstr(r*vars_f[s][i] <= expr);
            }   

            model.addConstr(vars_f[s][0] == 0);

            s++;
        }

        for (Vertex u = 0; u < num_vertices; u++) {
            for (int i = 0; i < B; i++) {
                GRBLinExpr expr = 0;

                for (int j = 0; j <= i; j++) {
                    expr += vars_x[u][j];
                }

                s = 0;
                for (auto comb : combs) {
                    auto depth = generate_covering_tree(g, r, comb);
                    
                    if (depth[u] - 1 < 0) {
                        s++;
                        continue; 
                    }

                    for (int j = 0; j <= i; j++) {
                        if (depth[u]-1 <= i-j) {
                            expr += vars_f[s][j];
                        }
                    }
                    s++;
                }
                model.addConstr(vars_y[u][i] <= expr);
            }
        }

        model.optimize();

        solve_status = inspect_gurobi_model(model);

        if (solve_status.status == GurobiSolveStatus::Optimal || 
            solve_status.status == GurobiSolveStatus::Feasible) {

            for(int i = 0; i < B; i++) {
                for(int u = 0; u < num_vertices; u++) {
                    if(vars_x[u][i].get(GRB_DoubleAttr_X) > 0.5) {
                        burning_sequence.push_back(u);
                    }
                }
            }
        } else {
            report_gurobi_failure(solve_status, model, "Cover");
        }
    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;

        report_gurobi_exception(solve_status, "Cover", e.getMessage());
    } catch (...) {
        cout << "Exception during optimization" << endl;
        report_gurobi_exception(solve_status, "Cover", "unknown exception during optimization");
    }
}

void find_burning_sequence(Graph const &g, int const r, int &LB, int &UB, vector<int> &burning_sequence, ggbp::ilp::GurobiParams params, GurobiStatus &solve_status) {
    PRECONDITION(r >= 1);
    PRECONDITION(g.order() >= r);

    vector<int> temp_burning_sequence;
    solve_status = {};

    int iteration = 1;
    while(LB <= UB) {
        int B = (UB + LB)/2;

        if (params.verbose) {
            cout << "Iteration #" << iteration << endl;
            cout << "Current lower bound = " << LB << endl;
            cout << "Current upper bound = " << UB << endl;
            cout << "Trying to find a burning sequence of length " << B << "." << endl;
        }

        params.time_limit -= solve_status.runtime_seconds;

        solve_ilp(g, r, B, params, temp_burning_sequence, solve_status);

        if (solve_status.error && solve_status.status != GurobiSolveStatus::Infeasible) {
            if (params.verbose)
                cout << "Error during ILP solving: " << solve_status.msg << endl;
            return;
        }


        if (solve_status.status == GurobiSolveStatus::Infeasible) {
            if (params.verbose) {
                cout << "There is no burning sequence of length " << B << "." << endl << endl;
            }
            LB = B+1;
            solve_status = {};
        } else {
            int turns = ggbp::utils::burning_turns(g, r, temp_burning_sequence);
            temp_burning_sequence.resize(turns);
            
            if (params.verbose) {
                cout << "A burning sequence of length " << turns << " was found." << endl << endl;
            }

            UB = turns;
            swap(burning_sequence, temp_burning_sequence);

            if (UB == LB)
                break;
        }

        iteration++;
    }
    
    if (params.verbose)
        cout << "Optimal burning sequence has length " << int(burning_sequence.size()) << "." << endl << endl;
}
}