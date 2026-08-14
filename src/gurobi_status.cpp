#include <string>
#include "gurobi_c++.h"

#include "ggbp/gurobi_status.hpp"

namespace ggbp::ilp {

GurobiSolveStatus classify_gurobi_status(
    const int native_status,
    const int solution_count) noexcept {
    if (native_status == GRB_OPTIMAL && solution_count > 0) {
        return GurobiSolveStatus::Optimal;
    }
    if (native_status == GRB_INFEASIBLE && solution_count == 0) {
        return GurobiSolveStatus::Infeasible;
    }
    if (solution_count > 0) {
        return GurobiSolveStatus::Feasible;
    }
    return GurobiSolveStatus::Unknown;
}


GurobiStatus inspect_gurobi_model(GRBModel& model) {
    const int native_status = model.get(GRB_IntAttr_Status);
    const int solution_count = model.get(GRB_IntAttr_SolCount);
    const double runtime_seconds = model.get(GRB_DoubleAttr_Runtime);
    std::optional<double> mip_gap;
    if (solution_count > 0 && model.get(GRB_IntAttr_IsMIP) != 0) {
        try {
            mip_gap = model.get(GRB_DoubleAttr_MIPGap);
        } catch (const GRBException&) {
            // Some specialized MIP solves do not expose a meaningful gap.
            // Missing diagnostic metadata must not discard a valid incumbent.
        }
    }
    return GurobiStatus{
        .status = classify_gurobi_status(native_status, solution_count),
        .native_status = native_status,
        .solution_count = solution_count,
        .runtime_seconds = runtime_seconds,
        .mip_gap = mip_gap,
    };
}

const char* gurobi_status_name(const int native_status) noexcept {
    switch (native_status) {
        case GRB_LOADED: return "LOADED";
        case GRB_OPTIMAL: return "OPTIMAL";
        case GRB_INFEASIBLE: return "INFEASIBLE";
        case GRB_INF_OR_UNBD: return "INF_OR_UNBD";
        case GRB_UNBOUNDED: return "UNBOUNDED";
        case GRB_CUTOFF: return "CUTOFF";
        case GRB_ITERATION_LIMIT: return "ITERATION_LIMIT";
        case GRB_NODE_LIMIT: return "NODE_LIMIT";
        case GRB_TIME_LIMIT: return "TIME_LIMIT";
        case GRB_SOLUTION_LIMIT: return "SOLUTION_LIMIT";
        case GRB_INTERRUPTED: return "INTERRUPTED";
        case GRB_NUMERIC: return "NUMERIC";
        case GRB_SUBOPTIMAL: return "SUBOPTIMAL";
        case GRB_INPROGRESS: return "INPROGRESS";
        case GRB_USER_OBJ_LIMIT: return "USER_OBJ_LIMIT";
        case GRB_WORK_LIMIT: return "WORK_LIMIT";
        case GRB_MEM_LIMIT: return "MEM_LIMIT";
        default: return "UNRECOGNIZED";
    }
}

void report_gurobi_failure(GurobiStatus& status, GRBModel& model, const std::string& context) {
    status.error = true;

    try {
        if (model.get(GRB_DoubleAttr_ObjVal) < status.ObjVal)
           status.ObjVal = model.get(GRB_DoubleAttr_ObjVal);
    } catch (...) { }

    try {
        status.ObjBound = model.get(GRB_DoubleAttr_ObjBound);
    } catch (...) { }
    
    status.msg =
        "Naive: Gurobi status " +
        std::string(gurobi_status_name(status.native_status)) +
        " without a usable incumbent";
}

void report_gurobi_exception(GurobiStatus& error, const std::string& context, const std::string& message) {
    error.error = true;
    error.status = GurobiSolveStatus::Error;
    error.native_status = 0;
    error.solution_count = 0;
    error.runtime_seconds = 0.0;
    error.mip_gap.reset();
    error.msg = context + ": " + message;
}
}