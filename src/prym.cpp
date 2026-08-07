#include <algorithm>
#include <exception>
#include <queue>
#include <stdexcept>
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

class callback_prym final : public GRBCallback {
public:
    callback_prym(
        const Graph& graph,
        const int horizon,
        const vector<vector<GRBVar>>& variables)
        : graph_(graph),
          n_vertices_(static_cast<int>(graph.order())),
          horizon_(horizon),
          vars_x_(variables) {}

    [[nodiscard]] bool failed() const noexcept {
        return callback_error_ != nullptr;
    }

    void rethrow_if_failed() const {
        if (callback_error_) {
            std::rethrow_exception(callback_error_);
        }
    }

protected:
    void callback() override {
        if (callback_error_) {
            abort();
            return;
        }

        try {
            if (where != GRB_CB_MIPSOL) {
                return;
            }

            const vector<int> burning_sequence = incumbent_sequence();
            add_violated_cover_constraint(burning_sequence);
        } catch (...) {
            callback_error_ = std::current_exception();
            abort();
        }
    }

private:
    Graph graph_;
    int n_vertices_;
    int horizon_;
    vector<vector<GRBVar>> vars_x_;
    std::exception_ptr callback_error_;

    vector<int> incumbent_sequence() {
        vector<int> sequence;
        sequence.reserve(static_cast<std::size_t>(horizon_));
        for (int turn = 0; turn < horizon_; ++turn) {
            for (int vertex = 0; vertex < n_vertices_; ++vertex) {
                if (getSolution(vars_x_[vertex][turn]) > 0.5) {
                    sequence.push_back(vertex);
                    break;
                }
            }
        }
        if (sequence.size() != static_cast<std::size_t>(horizon_)) {
            throw std::runtime_error(
                "PRYM callback received an incomplete incumbent sequence");
        }
        return sequence;
    }

    void add_violated_cover_constraint(const vector<int>& sequence) {
        for (const int source : sequence) {
            if (graph_.distance_matrix[source].empty()) {
                graph_.calculate_distances(source);
            }
        }

        int farthest_vertex = -1;
        int maximum_deficit = 0;
        for (int vertex = 0; vertex < n_vertices_; ++vertex) {
            int minimum_deficit = n_vertices_;
            for (int turn = 0; turn < horizon_; ++turn) {
                const int deficit =
                    graph_.get_distance(sequence[turn], vertex) -
                    (horizon_ - (turn + 1));
                minimum_deficit = std::min(minimum_deficit, deficit);
            }
            if (minimum_deficit > maximum_deficit) {
                maximum_deficit = minimum_deficit;
                farthest_vertex = vertex;
            }
        }

        if (farthest_vertex == -1) {
            return;
        }
        if (graph_.distance_matrix[farthest_vertex].empty()) {
            graph_.calculate_distances(farthest_vertex);
        }

        GRBLinExpr expression = 0;
        for (int turn = 0; turn < horizon_; ++turn) {
            for (int source = 0; source < n_vertices_; ++source) {
                if (graph_.get_distance(farthest_vertex, source) <=
                    horizon_ - (turn + 1)) {
                    expression += vars_x_[source][turn];
                }
            }
        }
        addLazy(expression >= 1);
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
        callback_prym cb(graph, B, vars_x);
        model.setCallback(&cb);
        model.optimize();
        cb.rethrow_if_failed();
        
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

    } catch (const GRBException& e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;

        report_gurobi_exception(err_info, "PRYM", e.getMessage());
        err_info.ObjVal = graph.order();
        err_info.ObjBound = 0;
    } catch (const std::exception& e) {
        report_gurobi_exception(err_info, "PRYM callback", e.what());
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
