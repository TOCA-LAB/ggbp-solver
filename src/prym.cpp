#include <queue>
#include <vector>
#include <cmath>
#include "gurobi_c++.h"
#include "ggbp/burning_simulation.hpp"
#include "ggbp/graph.hpp"
#include "ggbp/gurobi_status.hpp"

#include "ggbp/ggbp.hpp"

namespace ggbp::ilp::prym {
using namespace std;

void run_bff_d(Graph &graph, vector<int> &burning_sequence) {
    int n_vertices = graph.order();

    int max_vertex = 0, max_degree = graph.degree(0), curr_degree;

    for (int i = 1; i < n_vertices; i++) {
        curr_degree = graph.degree(i);
        if (curr_degree > max_degree) {
            max_vertex = i;
            max_degree = curr_degree;
        }
    }


    int n_burned = 0, selected, min_d, max_d, distance, u, to_burn;
    vector<int> queued(n_vertices, 0);
    queue<int> burning_queue;

    while (n_burned < n_vertices) {
        if (burning_sequence.size() == 0)
            selected = max_vertex;
        else {
            max_d = 0;
            for (int i = 0; i < n_vertices; i++) {
                min_d = n_vertices;
                for (int j : burning_sequence) {
                    distance = graph.get_distance(j, i);
                    if (distance < min_d)
                        min_d = distance;
                }
                if (min_d > max_d) {
                    max_d = min_d;
                    selected = i;
                }
            }
        }

        burning_sequence.push_back(selected);
        graph.calculate_distances(selected);
        
        if (!queued[selected]) {
            burning_queue.push(selected);
            queued[selected] = 1;
        }
        
        to_burn = burning_queue.size();
        while (to_burn > 0) {
            u = burning_queue.front();
            burning_queue.pop();
            n_burned++;

            for (int neighbor : graph.neighbors(u)) {
                if (!queued[neighbor]) {
                    queued[neighbor] = 1;
                    burning_queue.push(neighbor);
                }
            }
            to_burn--;
        }
    }
}

class callback_prym: public GRBCallback {
    public:
        Graph graph;
        int n_vertices;
        int B;
        vector<vector<GRBVar>> vars_x;
        callback_prym(Graph &_graph, int _B, vector<vector<GRBVar>> &_vars_x) {
            graph = _graph;
            n_vertices = graph.order();
            B = _B;
            vars_x = _vars_x;
        }
    protected:
        void callback() {
            try {
                if (where == GRB_CB_MIPSOL) {
                    vector<int> burning_sequence;
                    for (int j = 0; j < B; j++)
                        for (int i = 0; i < n_vertices; i++)
                            if (getSolution(vars_x[i][j]) > 0.5) {
                                burning_sequence.push_back(i);
                                break;   
                            }

                    for (int i : burning_sequence)
                        if (graph.distance_matrix[i].size() == 0)
                            graph.calculate_distances(i);

                    int min_d, max_v = -1, max_d = 0, diff = n_vertices;
                    for (int i = 0; i < n_vertices; i++) {
                        min_d = n_vertices;
                        for (int j = 0; j < B; j++) {
                            diff = graph.get_distance(burning_sequence[j], i) - (B - (j+1));
                            if (diff < min_d)
                                min_d = diff;
                        }
                        if (min_d > max_d) {
                            max_d = min_d;
                            max_v = i;
                        }
                    }

                    if (max_v != -1) {
                        if (graph.distance_matrix[max_v].size() == 0)
                            graph.calculate_distances(max_v);
                        GRBLinExpr expr = 0;
                        for (int j = 0; j < B; j++)
                            for (int u = 0; u < n_vertices; u++)
                                if (graph.get_distance(max_v, u) <= B - (j+1))
                                    expr += vars_x[u][j];
                        addLazy(expr >= 1);    
                    }

                }
            } catch (GRBException e) {
                cout << "Error number: " << e.getErrorCode() << endl;
                cout << e.getMessage() << endl;                
            } catch (...) {
                cout << "Error during callback" << endl;
            }
        }
};

void solve_prym(Graph &graph, int B, vector<int> &temp_burning_sequence, vector<int> initial_covered_vertices, GurobiParams params, GurobiErrorInfo &err_info) {
    temp_burning_sequence.clear();
    err_info = {};

    try {
        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, params.verbose);
        env.start();
        GRBModel model(env);

        model.set(GRB_DoubleParam_TimeLimit, params.time_limit);
        model.set(GRB_DoubleParam_SoftMemLimit, params.soft_mem_limit);
        model.set(GRB_DoubleParam_MemLimit, params.mem_limit);

        int n_vertices = graph.order();

        vector<vector<GRBVar>> vars_x(n_vertices, vector<GRBVar>(B));
        for (int i = 0; i < n_vertices; i++)
            for (int j = 0; j < B; j++)
                vars_x[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "x_" + to_string(i) + "_" + to_string(j));

        for (int i = 0; i < n_vertices; i++) {
            GRBLinExpr expr = 0;
            for (int j = 0; j < B; j++)
                expr += vars_x[i][j];
            model.addConstr(expr <= 1);
        }

        for (int j = 0; j < B; j++) {
            GRBLinExpr expr = 0;
            for (int i = 0; i < n_vertices; i++)
                expr += vars_x[i][j];
            model.addConstr(expr == 1);
        }

        for (int i : initial_covered_vertices) {
            GRBLinExpr expr = 0;
            for (int j = 0; j < B; j++)
                for (int u = 0; u < n_vertices; u++)
                    if (graph.get_distance(i, u) <= B - (j+1))
                        expr += vars_x[u][j];
            model.addConstr(expr >= 1);
        }

        for(int j = 0; j < B; j++)
            for(int i = 0; i < n_vertices; i++)
                vars_x[i][j].set(GRB_IntAttr_BranchPriority, B-j);
                
        model.set(GRB_IntParam_LazyConstraints, 1);
        callback_prym cb = callback_prym(graph, B, vars_x);
        model.setCallback(&cb);
        model.optimize();
        
        const auto solve_status = inspect_gurobi_model(model);
        record_gurobi_status(err_info, solve_status);
        if (solve_status.status == GurobiSolveStatus::Optimal ||
            solve_status.status == GurobiSolveStatus::Feasible) {
            for(int j = 0; j < B; j++) {
                for(int i = 0; i < n_vertices; i++) {
                    if(vars_x[i][j].get(GRB_DoubleAttr_X) > 0.5) {
                        temp_burning_sequence.push_back(i);
                    }
                }
            }

            const auto certificate = simulation::simulate_burning(
                graph, 1, temp_burning_sequence);
            if (!certificate.feasible || !certificate.completion_turns) {
                temp_burning_sequence.clear();
                err_info.error = true;
                err_info.status = GurobiSolveStatus::Error;
                err_info.msg =
                    "PRYM returned an invalid burning-sequence incumbent";
                return;
            }
            temp_burning_sequence.resize(*certificate.completion_turns);
            if (solve_status.status == GurobiSolveStatus::Feasible) {
                err_info.msg = "PRYM accepted a feasible incumbent at " +
                    string(gurobi_status_name(solve_status.native_status));
            }
        } else {
            report_gurobi_failure(err_info, solve_status, "PRYM");
        }

    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;

        report_gurobi_exception(err_info, "PRYM", e.getMessage());
        err_info.ObjVal = graph.order();
        err_info.ObjBound = 0;
    } catch (...) {
        cout << "Exception during optimization" << endl;
        report_gurobi_exception(
            err_info, "PRYM", "unknown exception during optimization");
        err_info.ObjVal = graph.order();
        err_info.ObjBound = 0;
    }
}


void run_prym(Graph &graph, vector<int> &burning_sequence, GurobiParams params, GurobiErrorInfo &err_info) {
    int LB, UB, B;
    vector<int> temp_burning_sequence, initial_covered_vertices;

    run_bff_d(graph, temp_burning_sequence);

    burning_sequence = temp_burning_sequence;
    initial_covered_vertices = temp_burning_sequence;

    UB = int(temp_burning_sequence.size());
    LB = (UB + 4) / 3; 


    if (params.verbose)
        cout << "BFF-d found a burning sequence of length " << UB << endl << endl;
    
    while (LB < UB) {
        B = (LB + UB)/2;

        if (params.verbose) {
            cout << "Current lower bound = " << LB << endl;
            cout << "Current upper bound = " << UB << endl;
            cout << "Trying to find a burning sequence of length " << B << "." << endl << endl;
        }

        solve_prym(graph, B, temp_burning_sequence, initial_covered_vertices, params, err_info);

        if (err_info.status == GurobiSolveStatus::Infeasible) {
            err_info = {};
            temp_burning_sequence.clear();
        } else if (err_info.error) {
            return;
        }

        if (params.verbose)
        cout << temp_burning_sequence.size() << endl;

        if (temp_burning_sequence.size() > 0) {
            swap(burning_sequence, temp_burning_sequence);
            UB = B;
        } else {
            LB = B + 1;
        }
    }
}

}
