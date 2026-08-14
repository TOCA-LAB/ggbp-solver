#ifndef GUROBI_STATUS_H
#define GUROBI_STATUS_H

#include <optional>
#include <string>
#include "gurobi_c++.h"

namespace ggbp::ilp {
    enum class GurobiSolveStatus {
        NotStarted,
        Optimal,
        Feasible,
        Infeasible,
        Unknown,
        Error,
    };

    struct GurobiStatus {
        bool error = false;
        GurobiSolveStatus status = GurobiSolveStatus::NotStarted;
        int native_status = 0;
        int solution_count = 0;
        double runtime_seconds = 0.0;
        std::optional<double> mip_gap;
        double ObjVal = 0.0;
        double ObjBound = 0.0;
        std::string msg;
    };

    GurobiStatus inspect_gurobi_model(GRBModel& model);

    [[nodiscard]] GurobiSolveStatus classify_gurobi_status(int native_status, int solution_count) noexcept;
    [[nodiscard]] const char* gurobi_status_name(int native_status) noexcept;

    void report_gurobi_failure(GurobiStatus& status, GRBModel& model, const std::string& context);

    void report_gurobi_exception(GurobiStatus& error, const std::string& context, const std::string& message);
}
#endif
