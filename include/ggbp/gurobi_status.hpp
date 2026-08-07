#ifndef GGBP_GUROBI_STATUS_HPP
#define GGBP_GUROBI_STATUS_HPP

#include <optional>
#include <string>

#include <gurobi_c++.h>

#include "ggbp/ggbp.hpp"

namespace ggbp::ilp {

struct GurobiStatusSnapshot {
    GurobiSolveStatus status = GurobiSolveStatus::NotStarted;
    int native_status = 0;
    int solution_count = 0;
    double runtime_seconds = 0.0;
    std::optional<double> mip_gap;

    [[nodiscard]] bool has_incumbent() const noexcept {
        return solution_count > 0;
    }
};

[[nodiscard]] GurobiSolveStatus classify_gurobi_status(
    int native_status,
    int solution_count) noexcept;

[[nodiscard]] GurobiStatusSnapshot inspect_gurobi_model(GRBModel& model);

[[nodiscard]] const char* gurobi_status_name(int native_status) noexcept;

void record_gurobi_status(
    GurobiErrorInfo& error,
    const GurobiStatusSnapshot& snapshot);

void report_gurobi_failure(
    GurobiErrorInfo& error,
    const GurobiStatusSnapshot& snapshot,
    const std::string& context);

void report_gurobi_exception(
    GurobiErrorInfo& error,
    const std::string& context,
    const std::string& message);

} // namespace ggbp::ilp

#endif
