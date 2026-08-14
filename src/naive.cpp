#include <cmath>
#include <optional>
#include "gurobi_c++.h"
#include "utils/contract.hpp"
#include "ggbp/graph.hpp"
#include "utils/utils.hpp"
#include "ggbp/gurobi_status.hpp"

#include "ggbp/ggbp.hpp"

namespace ggbp::ilp {
using namespace std;

int get_lower_bound(Graph const &g, int const r, bool const verbose) {
    int count = 0;
    auto const threshold = r <= 0 ? Graph::size_type{0} : static_cast<Graph::size_type>(r);
    
    for (Graph::size_type u = 0; u < g.order(); u++) {
        if (g.degree(u) < threshold)
            count++;
    }

    return count;
}

int get_upper_bound(
    Graph const &graph,
    int const r,
    bool const verbose,
    unsigned const seed) {
    ggbp::brkga::BRKGAParams params{};
    
    params.p = 60;
    params.MAX_GENS = 25;
    params.pe = 0.20;
    params.pm = 0.10;
    params.rhoe = 0.70;	
    params.K = 2;
    params.MAXT = 2;
    params.seed = seed;
    params.verbose = verbose;

    vector<int> burning_sequence, iters;
    int first_best;

    ggbp::brkga::brkga(graph, r, params, burning_sequence, iters, first_best);


    return burning_sequence.size();
}
}

namespace ggbp::ilp::naive {
using namespace std;

void add_radius_cut(GRBModel &model, Graph &g, vector<vector<GRBVar>> x, vector<vector<GRBVar>> y, int const r, int const B) {
    auto num_vertices = g.order();

    for (int u = 0; u < num_vertices; u++) {
        for (int t = 1; t < B+1; t++) {
            GRBLinExpr expr = 0;
            for (int v = 0; v < num_vertices; v++) {
                for (int ti = 1; ti <= t; ti++) {
                    if (g.get_distance(u, v) <= (t - ti))
                        expr += y[v][ti];
                }
            }
            model.addConstr(x[u][t] <= expr, "Radius_" + to_string(u) + "_" + to_string(t));
        }
    }
}


void add_active_neighborhood_logic(GRBModel &model, Graph const &g, vector<vector<GRBVar>> x, vector<vector<GRBVar>> w, int const r, int const B) {
    auto num_vertices = g.order();

    for (int u = 0; u < num_vertices; u++) {
        if (g.degree(u) < r)
            continue;

        for (int t = 2; t < B+1; t++) {
            GRBLinExpr expr = 0;

            for (auto v : g.neighbors(u)) {
                expr += x[v][t-1] - x[v][t-2];
            }

            model.addConstr(w[u][t] <= expr, "Ngbr_" + to_string(u) + "_" + to_string(t));
        }
    }
}

std::optional<double> solve_ilp(Graph &g, int const r, int const B, LPParams params, std::vector<int> &burning_sequence, GurobiStatus &solve_status) {
    int const n_vertices = g.order();

    burning_sequence.clear();

    solve_status = {};
    solve_status.ObjVal = B;
    solve_status.ObjBound = 0;

    char vtype = GRB_BINARY;

    if (params.linear_relax) {
        vtype = GRB_CONTINUOUS;
    }

    try {        
        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, params.verbose);
        env.start();
        GRBModel model(env);

        model.set(GRB_DoubleParam_TimeLimit, params.time_limit);
        model.set(GRB_DoubleParam_SoftMemLimit, params.soft_mem_limit);
        model.set(GRB_DoubleParam_MemLimit, params.mem_limit);

        // --- Variáveis ---
        // x[v, t]: 1 se o vértice v está queimado no tempo t
        vector<vector<GRBVar>> x(n_vertices, vector<GRBVar>(B+1));

        // y[v, t]: 1 se o vértice v é inflamado manualmente no tempo t
        vector<vector<GRBVar>> y(n_vertices, vector<GRBVar>(B+1));
        
        // w[v, t]: 1 se o vértice v é inflamado espontaneamente no tempo t
        vector<vector<GRBVar>> w(n_vertices, vector<GRBVar>(B+1));

        for (int u = 0; u < n_vertices; u++) {
            x[u][0] = model.addVar(0, 1, 0, vtype, "x_" + to_string(u) + "_0");

            for (int t = 1; t < B + 1; t++) {
                x[u][t] = model.addVar(0, 1, 0, vtype, "x_" + to_string(u) + "_" + to_string(t));
                y[u][t] = model.addVar(0, 1, 0, vtype, "y_" + to_string(u) + "_" + to_string(t));
                w[u][t] = model.addVar(0, 1, 0, vtype, "w_" + to_string(u) + "_" + to_string(t));
            }
        }

        // --- Restrições ---
        // 0. Estado inicial: ninguém queimado
        for (int u = 0; u < n_vertices; u++) {
            model.addConstr(x[u][0] == 0, "init");
        }

        GRBLinExpr expr;
        for (int t = 1; t < B + 1; t++) {
            // 1. Exatamente uma ignição manual por passo
            if (!params.opt) {
                expr = 0;
                for (int u = 0; u < n_vertices; u++) {
                    expr += y[u][t];
                }
                model.addConstr(expr == 1, "ignitor_" + to_string(t));
            }

            for (int u = 0; u < n_vertices; u++) {
                // 2. Persistência: uma vez queimado, sempre queimado
                model.addConstr(x[u][t] >= x[u][t-1], "stay_" + to_string(u) + "_" + to_string(t));
            
                //  3. Lógica de Difusão (Threshold r)
                if (g.degree(u) >= r) {
                    // w[v,t] só pode ser 1 se r vizinhos estavam queimados em t-1
                    expr = 0;
                    for (int v : g.neighbors(u))
                        expr += x[v][t-1];
                    model.addConstr(r * w[u][t] <= expr, "diff_" + to_string(u) + "_" + to_string(t));
                } else {
                    model.addConstr(w[u][t] == 0);
                }

                // 4. Transição: x[v,t] é 1 se já estava queimado, foi ignitado ou difundido
                model.addConstr(x[u][t] <= x[u][t-1] + y[u][t] + w[u][t], "trans_" + to_string(u) + "_" + to_string(t));
            }
        }

        // 5. Condição de Finalização: TODOS os vértices devem estar queimados no tempo B
        for (int u = 0; u < n_vertices; u++) {
            model.addConstr(x[u][B] == 1, "burned_" + to_string(u));
        }

        if (params.active_neighborhood)
            add_active_neighborhood_logic(model, g, x, w, r, B);

        if (params.radius_cut)
            add_radius_cut(model, g, x, y, r, B);

        // vector<Vertex> sequence = {2, 33, 5, 6, 11};

        // for (int t = 1; t < B; t++) {
        //     model.addConstr(y[sequence[t-1]][t] == 1);
        // }

        // Se escolher versão de otimização
        if (params.opt) {
            // Variável de otimização
            // done[t]: 1 se todos os vértices estão queimados no turno t
            vector<GRBVar> done(B+1);

            for (int t = 0; t < B+1; t++) {
                done[t] = model.addVar(0, 1, 0, vtype, "done_" + to_string(t));
            }

            for (int t = 1; t < B+1; t++) {
                expr = 0;
                for (int u = 0; u < n_vertices; u++) {
                    expr += y[u][t];
                }
                model.addConstr(expr == 1 - done[t-1], "ignition_limit_done_" + to_string(t));

                model.addConstr(done[t] >= done[t-1], "done_persistence" + to_string(t));

                for (int u = 0; u < n_vertices; u++) {
                    model.addConstr(done[t] <= x[u][t], "check_done_" + to_string(u) + "_" + to_string(t));
                }
            }
            
            // Minimiza a quantidade de turnos até todos vértices serem queimados
            GRBLinExpr obj = 0;
            for (int t = 1; t < B+1; t++) {
                obj += 1 - done[t];
            }
            model.setObjective(obj + 1, GRB_MINIMIZE);
        }


        model.optimize();

        solve_status = inspect_gurobi_model(model);

        if (solve_status.status == GurobiSolveStatus::Optimal || 
            solve_status.status == GurobiSolveStatus::Feasible) {
            
            if (params.linear_relax) {
                if (solve_status.status == GurobiSolveStatus::Optimal)
                    return model.get(GRB_DoubleAttr_ObjVal);

                solve_status.error = true;
                solve_status.ObjVal = model.get(GRB_DoubleAttr_ObjVal);
                solve_status.ObjBound = model.get(GRB_DoubleAttr_ObjBound);
                solve_status.msg =
                    "Naive LP stopped with a non-optimal feasible point at " +
                    string(gurobi_status_name(solve_status.native_status));
                return std::nullopt;
            }

            for(int t = 1; t < B+1; t++) {
                for(int u = 0; u < n_vertices; u++) {
                    if (y[u][t].get(GRB_DoubleAttr_X) > 0.5) {
                        burning_sequence.push_back(u);
                    }
                }
            }
            
            if (solve_status.status == GurobiSolveStatus::Feasible &&
                params.opt) {
                solve_status.error = true;
                solve_status.ObjVal = model.get(GRB_DoubleAttr_ObjVal);
                solve_status.ObjBound = model.get(GRB_DoubleAttr_ObjBound);
                solve_status.msg =
                    "Naive retained a non-optimal integer incumbent at " +
                    string(gurobi_status_name(solve_status.native_status));
            } else if (solve_status.status == GurobiSolveStatus::Feasible) {
                solve_status.msg = "Naive accepted a feasible incumbent at " +
                    string(gurobi_status_name(solve_status.native_status));
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

    return std::nullopt;
}

void find_burning_sequence(Graph &g, int const r, int &LB, int &UB, LPParams params, std::vector<int> &burning_sequence, GurobiStatus &solve_status) {
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
