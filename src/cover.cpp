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

void solve_ilp(Graph const &g, int const r, int const B, vector<int> &burning_sequence, ggbp::ilp::GurobiParams params, GurobiErrorInfo &err_info) {
    int num_vertices = g.order();
    vector<int> vertices = g.vertices();

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

        // for (auto u = 0; u < g.order(); u++) {
        //     if (g.degree(u) < r) {
        //         GRBLinExpr expr = 0;
                
        //         for (int i = 0; i < B; i++) {
        //             expr += vars_x[u][i];
        //         }
        //         model.addConstr(expr == 1);
        //         cout << u << endl;
        //     }
        // }

        // vector<int> seq =  {1, 3, 4};
        // for (int i = 0; i < seq.size(); i++) {
        //     model.addConstr(vars_x[seq[i]][i] == 1);
        // }
    
        // GRBLinExpr obj = 0;
        // for (int u = 0; u < num_vertices; u++) {
        //     obj += vars_y[u][B-1];
        // }
        // model.setObjective(obj, GRB_MAXIMIZE);

        model.optimize();

        int opt_status = model.get(GRB_IntAttr_Status);
        burning_sequence.clear();

        if (opt_status == GRB_OPTIMAL) {
            for(int i = 0; i < B; i++) {
                for(int u = 0; u < num_vertices; u++) {
                    if(vars_x[u][i].get(GRB_DoubleAttr_X) > 0.5) {
                        burning_sequence.push_back(u);
                    }
                    
                }
            }
            Vertex u;
            cout << ggbp::utils::is_burning_sequence(g, r, burning_sequence, u) << endl;

            for (auto v : burning_sequence) {
                cout << v << " ";
            }
            cout << endl;
        } else if (opt_status == GRB_TIME_LIMIT) {
            err_info.error = true;
            
            if (model.get(GRB_DoubleAttr_ObjVal) < B)
                err_info.ObjVal = model.get(GRB_DoubleAttr_ObjVal);
            else
                err_info.ObjVal = B;
            
            err_info.ObjBound = model.get(GRB_DoubleAttr_ObjBound);
            err_info.msg = "time limit";
        } else if (opt_status == GRB_MEM_LIMIT) {
            err_info.error = true;

            if (model.get(GRB_DoubleAttr_ObjVal) < B)
                err_info.ObjVal = model.get(GRB_DoubleAttr_ObjVal);
            else
                err_info.ObjVal = B;

            err_info.ObjBound = model.get(GRB_DoubleAttr_ObjBound);
            err_info.msg = "memory limit";
        }
    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;

        if (e.getErrorCode() == GRB_ERROR_OUT_OF_MEMORY) {
            err_info.error = true;
            err_info.ObjVal = B;
            err_info.ObjBound = 0;
            err_info.msg = "out of memory";
        }
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }
}

void find_burning_sequence(Graph const &g, int const r, int &LB, int &UB, vector<int> &burning_sequence, ggbp::ilp::GurobiParams params, GurobiErrorInfo &err_info) {
    PRECONDITION(r >= 1);
    PRECONDITION(G.order() >= r);

    vector<int> temp_burning_sequence;

    int iteration = 1;
    while(LB <= UB) {
        int B = (UB + LB)/2;

        if (params.verbose) {
            cout << "Iteration #" << iteration << endl;
            cout << "Current lower bound = " << LB << endl;
            cout << "Current upper bound = " << UB << endl;
            cout << "Trying to find a burning sequence of length " << B << "." << endl;
        }

        solve_ilp(g, r, B, temp_burning_sequence, params, err_info);

        if (err_info.error)
            return;

        int turns = ggbp::utils::burning_turns(g, r, temp_burning_sequence);
        // temp_burning_sequence.resize(turns);

        if(int(temp_burning_sequence.size()) == 0) {
            if (params.verbose) {
                cout << "There is no burning sequence of length " << B << "." << endl << endl;
            }
            LB = B+1;
        } else {
            if (params.verbose) {
                cout << "A burning sequence of length " << B << " was found." << endl << endl;
            }
            UB = B;
            swap(burning_sequence, temp_burning_sequence);

            if (UB == LB)
                break;
        }

        iteration++;
    }

    cout << "Optimal burning sequence has length " << int(burning_sequence.size()) << "." << endl << endl;
}
}