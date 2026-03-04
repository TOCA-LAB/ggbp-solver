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
            tree_neighbors[v]++;

            if (tree_neighbors[v] >= r && depth[v] == -1) {
                depth[v] = depth[u] + 1;
                q.push(v);
            }
            // g.remove_edge(u, v);
        }
    }

    return depth;
}


class LazyCovering: public GRBCallback {
    public:
    Graph const &g;
    vector<vector<GRBVar>> &vars_x;
    vector<vector<GRBVar>> &vars_y;
    vector<vector<GRBVar>> &vars_f;
    vector<bool> &covered;
    int B;
    int r;

    LazyCovering(Graph const &g, int r, int B, vector<vector<GRBVar>> &vars_x,
                 vector<vector<GRBVar>> &vars_y, vector<vector<GRBVar>> &vars_f,
                 vector<bool> &covered)
    : g(g), vars_x(vars_x), vars_y(vars_y), vars_f(vars_f), covered(covered)
    {
        this->r = r;
        this->B = B;
    }

    protected:
    void callback() {
        try {
            if (where == GRB_CB_MIPSOL) {
                cout << "callback: mipsol" << endl;
                int num_vertices = g.order();
                vector<int> burning_sequence;

                for(int i = 0; i < B; i++) {
                    for(Vertex u = 0; u < num_vertices; u++) {
                        if (getSolution(vars_x[u][i]) > 0.5) {
                            burning_sequence.push_back(u);
                            u = num_vertices;
                        }
                    }
                }

                for (auto u : burning_sequence)
                    cout << u << " ";
                cout << endl;

                Vertex unburned;
                if (!ggbp::utils::is_burning_sequence(g, r, burning_sequence, unburned)) {
                    cout << "Lazy vertex: " << unburned << endl;

                    queue<Vertex> q;
                    q.push(unburned); 

                    auto vertices = g.vertices();
                    auto combs = ggbp::utils::Combinations(vertices, r);
                    
                    while (!q.empty()) {
                        Vertex u = q.front();
                        q.pop();

                        for (int i = 0; i < B; i++) {
                            GRBLinExpr expr = 0;
                            double total = 0;

                            for (int j = 0; j <= i; j++) {
                                expr += vars_x[u][j];
                                total += getSolution(vars_x[u][i]);
                            }

                            unsigned long long s = 0;
                            for (auto comb : combs) {
                                auto depth = generate_covering_tree(g, r, comb);
                                
                                if (depth[u] - 1 < 0) {
                                    s++;
                                    continue; 
                                }

                                for (int j = 0; j <= i; j++) {
                                    if (depth[u] - 1 <= i-j) {
                                        expr += vars_f[s][j];

                                        if (getSolution(vars_f[s][i]) > 0.5) {
                                            cout << "(";
                                
                                            for (auto c : comb)
                                                cout << c << ", ";

                                            cout << ")/" << i << " : d[u] = " << depth[u] << " / var = "
                                                << getSolution(vars_f[s][i]) << endl;
                                        }
                                        total += getSolution(vars_f[s][i]);
                                    }
                                }
                                s++;
                            }
                            cout << total << " / " << (1 <= total) << endl;
                            addLazy(vars_y[unburned][i] <= expr);
                        }
                    }
                    // unsigned long long s = 0;
                    // for (auto comb : combs) {
                    //     auto depth = generate_covering_tree(g, r, comb);

                    //     if (depth[unburned] == -1) {
                    //         s++;
                    //         continue; 
                    //     }

                    //     for (int i = 0; i < B; i++) {
                    //         if (depth[unburned] < (B-i)) {
                    //             expr += vars_f[s][i];  
                    //             if (getSolution(vars_f[s][i]) > 0.5) {
                    //                 cout << "(";
                        
                    //                 for (auto c : comb)
                    //                     cout << c << ", ";

                    //                 cout << ")/" << i << " : d[u] = " << depth[unburned] << " / var = "
                    //                      << getSolution(vars_f[s][i]) << endl;
                    //             }
                    //             total += getSolution(vars_f[s][i]);
                    //         }
                    //     }
                    //     s++;
                    // }
                    // cout << total << " / " << (1 <= total) << endl;
                    // addLazy(1 <= expr);
                    cout << "Added Lazy" << endl << endl;
                }

            }
            // cout << endl;
        } catch (GRBException e) {
            cout << "Error number: " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch (...) {
            cout << "Error during callback" << endl;	
        }
    }
};

void solve_covering_ilp(Graph const &g, int const r, int const B, vector<int> &burning_sequence, bool verbose) {
    int num_vertices = g.order();
    vector<int> vertices = g.vertices();

    try {
        if (verbose)
            cout << endl;
        
        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, verbose);
        env.start();

        GRBModel model(env);
        
        model.set(GRB_DoubleParam_TimeLimit, 7200.0);
        model.set(GRB_DoubleParam_SoftMemLimit, 13.0);
        model.set(GRB_DoubleParam_MemLimit, 13.5);

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
                model.addConstr(vars_y[u][i] <= expr).set(GRB_IntAttr_Lazy, 2);
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

        vector<bool> covered(num_vertices, false);

        // model.set(GRB_IntParam_LazyConstraints, 1);
        // LazyCovering cb = LazyCovering(g, r, B, vars_x, vars_y, vars_f, covered);
        // model.setCallback(&cb);
        model.optimize();

        int opt_status = model.get(GRB_IntAttr_Status);
        burning_sequence.clear();
        if (opt_status == GRB_OPTIMAL) {
            for(int i = 0; i < B; i++) {
                for(int u = 0; u < num_vertices; u++) {
                    if(vars_x[u][i].get(GRB_DoubleAttr_X) > 0.5) {
                        burning_sequence.push_back(u);
                        // u = num_vertices;
                    }
                    
                    // cout << vars_x[u][i].get(GRB_StringAttr_VarName) << ": " << (1 ? vars_x[u][i].get(GRB_DoubleAttr_X) > 0.5 : 0) << endl;
                    // cout << vars_y[u][i].get(GRB_StringAttr_VarName) << ": " << (1 ? vars_y[u][i].get(GRB_DoubleAttr_X) > 0.5 : 0);
                    
                    // double sum = 0;

                    // for (int j = 0; j <= i; j++) {
                    //     sum += vars_x[u][j].get(GRB_DoubleAttr_X);
                    // }

                    // s = 0;
                    // for (auto comb : combs) {
                    //     auto depth = generate_covering_tree(g, r, comb);
                        
                    //     if (depth[u] - 1 >= 0) {
                    //         s++;
                    //         continue; 
                    //     }

                    //     for (int j = 0; j < i; j++) {
                    //         if (depth[u]-1 < i-j) {
                    //             sum += vars_f[s][j].get(GRB_DoubleAttr_X);
                    //             if (vars_f[s][j].get(GRB_DoubleAttr_X) > 0.5) {
                    //                 cout << "(";
                    //                 for (auto v : comb)
                    //                     cout << v << " ";
                    //                 cout << ")";
                    //             }
                    //         }
                    //     }
                    //     s++;
                    // }
                    // cout << " " << sum << endl;
                    
                }
            }

            if (!burning_sequence.empty()) {
                Vertex unburned;
                if (!ggbp::utils::is_burning_sequence(g, r, burning_sequence, unburned))
                    cout << "Solution is not burning sequence" << endl;
            }

            // unsigned long long s = 0;
            // for (auto comb : combs) {
            //     for (int i = 0; i < B; i++) {
            //         if (vars_f[s][i].get(GRB_DoubleAttr_X) > 0.5) {
            //             cout << vars_f[s][i].get(GRB_StringAttr_VarName) << ": 1 (";

            //             for (auto u : comb)
            //                 cout << u << " ";
            //             cout << ")" << endl;
            //         }
            //         double sum = 0;
            //         for (auto u : s) {
            //             sum += vars_y[u][j].get(GRB_DoubleAttr_X);
            //         cout << " " << sum << endl;
            //         }
            //     }
            //     s++;
            // }
        }
    } catch (GRBException e) {
        
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }
}

int get_lower_bound(Graph const &g, int const r) {
    int count = 0;
    
    for (auto u = 0; u < g.order(); u++) {
        if (g.degree(u) < r)
            count++;
    }

    return 0;
}

int get_upper_bound(Graph const &g, int const r) {
    // return g.order();
    return 15;
}

void ggbp_covering(Graph const &g, int const r, vector<int> &burning_sequence, bool verbose) {
    PRECONDITION(r >= 1);
    PRECONDITION(G.order() >= r);

    int LB = get_lower_bound(g, r);
    int UB = get_upper_bound(g, r);
    vector<int> temp_burning_sequence;

    int iteration = 0;
    while(LB <= UB) {
        int B = (UB + LB)/2;

        if (1) {
            cout << "Iteration #" << iteration << endl;
            cout << "Current lower bound = " << LB << endl;
            cout << "Current upper bound = " << UB << endl;
            cout << "Trying to find a burning sequence of length " << B << "." << endl;
        }

        solve_covering_ilp(g, r, B, temp_burning_sequence, verbose);

        if(int(temp_burning_sequence.size()) == 0) {
            if (1) {
                cout << "There is no burning sequence of length " << B << "." << endl << endl;
            }
            LB = B+1;
        } else {
            if (1) {
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