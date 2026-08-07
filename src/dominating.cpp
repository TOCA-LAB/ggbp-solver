#include <vector>
#include <optional>
#include <string>
#include <iostream>
#include <exception>
#include <stdexcept>
#include "ggbp/burning_simulation.hpp"
#include "ggbp/graph.hpp"
#include "gurobi_c++.h"

#include "ggbp/ggbp.hpp"

namespace ggbp::heuristics {
using std::vector;

void min_dominating_set(Graph const &graph, vector<int> &dominating_set, ggbp::ilp::GurobiParams params, ggbp::ilp::GurobiErrorInfo &err_info) {
    try {
        int n_vertices = graph.order();

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
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;

        err_info.error = true;
        err_info.ObjVal = graph.order();
        err_info.ObjBound = 0;
        if (e.getErrorCode() == GRB_ERROR_OUT_OF_MEMORY) {
            err_info.msg = "out of memory";
        }
    } catch (...) {
        std::cout << "Exception during optimization" << std::endl;
        err_info.error = true;
        err_info.ObjVal = graph.order();
        err_info.ObjBound = 0;
    }
}

namespace {

using simulation::BurningState;

std::optional<Vertex> next_preferred_source(
    const vector<Vertex>& preferred_sources,
    std::size_t& next_preferred,
    const BurningState& state) {
    while (next_preferred < preferred_sources.size()) {
        const Vertex source = preferred_sources[next_preferred++];
        if (!state.is_burned(source) && !state.is_in_frontier(source)) {
            return source;
        }
    }
    return std::nullopt;
}

std::size_t newly_scheduled_neighbors(
    const Graph& graph,
    const Vertex source,
    const BurningState& state) {
    std::size_t count = 0;
    for (const Vertex neighbor : graph.neighbors(source)) {
        if (!state.is_burned(neighbor) &&
            !state.is_in_frontier(neighbor) &&
            state.burned_neighbor_count(neighbor) + 1 ==
                static_cast<unsigned>(state.threshold())) {
            ++count;
        }
    }
    return count;
}

std::optional<Vertex> best_completion_source(
    const Graph& graph,
    const BurningState& state) {
    std::optional<Vertex> best_source;
    std::size_t best_score = 0;

    for (Vertex vertex = 0;
         static_cast<Graph::size_type>(vertex) < graph.order();
         ++vertex) {
        if (state.is_burned(vertex) || state.is_in_frontier(vertex)) {
            continue;
        }

        const std::size_t score =
            newly_scheduled_neighbors(graph, vertex, state);
        if (!best_source || score > best_score) {
            best_source = vertex;
            best_score = score;
        }
    }
    return best_source;
}

std::optional<Vertex> first_frontier_vertex(const BurningState& state) {
    for (const Vertex vertex : state.frontier()) {
        if (!state.is_burned(vertex)) {
            return vertex;
        }
    }
    return std::nullopt;
}

} // namespace

vector<Vertex> dominating(
    const Graph& graph,
    const int r,
    const vector<Vertex>& dominating_set) {

    BurningState state(graph, r);
    for (const Vertex vertex : dominating_set) {
        (void)state.is_burned(vertex);
    }

    vector<Vertex> burning_sequence;
    burning_sequence.reserve(graph.order());
    std::size_t next_preferred = 0;

    while (!state.complete()) {
        auto source = next_preferred_source(
            dominating_set, next_preferred, state);
        if (!source) {
            source = best_completion_source(graph, state);
        }
        if (!source) {
            source = first_frontier_vertex(state);
        }
        if (!source) {
            throw std::logic_error(
                "Dominating could not select an unburned source");
        }

        burning_sequence.push_back(*source);
        state.burn(*source);
    }

    if (!simulation::simulate_burning(graph, r, burning_sequence).feasible) {
        throw std::logic_error(
            "Dominating produced an invalid burning sequence");
    }
    return burning_sequence;
}

vector<Vertex> dominating(
    const Graph& graph,
    const int r,
    const ggbp::ilp::GurobiParams params,
    ggbp::ilp::GurobiErrorInfo& err_info) {
    vector<int> dominating_set;
    err_info = {};
    min_dominating_set(graph, dominating_set, params, err_info);

    if (err_info.error) {
        return {};
    }
    return dominating(graph, r, dominating_set);
}
}
