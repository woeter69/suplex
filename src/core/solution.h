#pragma once

#include "types.h"
#include <vector>
#include <cstdint>

namespace suplex {

/**
 * @brief Container for solver outputs, including primal/dual variables,
 * objective values, and runtime statistics.
 */
struct Solution {
    SolverStatus status = SolverStatus::NOT_STARTED;
    Real objective_value = 0.0;
    std::vector<Real> primal_values;     // x (size: num_cols)
    std::vector<Real> dual_values;       // y (size: num_rows)
    std::vector<Real> reduced_costs;     // rc (size: num_cols)
    std::vector<Real> row_activities;    // Ax (size: num_rows)

    // MILP-specific fields
    Real best_bound = 0.0;
    Real mip_gap = 0.0;
    int64_t num_nodes = 0;

    // Execution statistics
    int64_t num_iterations = 0;
    double solve_time_seconds = 0.0;

    Solution() = default;
    explicit Solution(SolverStatus s) : status(s) {}

    bool is_optimal() const { return status == SolverStatus::OPTIMAL; }
    bool is_feasible() const;
    void reset();
};

} // namespace suplex
