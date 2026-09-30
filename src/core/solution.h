#pragma once
#include "types.h"
#include <vector>
#include <cstdint>

namespace suplex {

// ── Solution container ────────────────────────────────────────────────────────
///
/// Holds the complete result of a solve call.  Always check `status` first.
///
struct Solution {
    SolverStatus status = SolverStatus::NOT_STARTED;

    // Primal / dual
    Real              objective_value = 0.0;
    std::vector<Real> primal_values;    ///< x  — size num_cols
    std::vector<Real> dual_values;      ///< y  — size num_rows (shadow prices)
    std::vector<Real> reduced_costs;    ///< rc — size num_cols
    std::vector<Real> row_activities;   ///< Ax — size num_rows

    // MILP-specific
    Real     best_bound = 0.0;   ///< Best dual bound from open B&B nodes
    Real     mip_gap    = 1.0;   ///< Relative MIP gap: |obj - bound| / |obj|
    int64_t  num_nodes  = 0;     ///< B&B nodes explored

    // Statistics
    int64_t  num_iterations    = 0;
    double   solve_time_seconds = 0.0;

    // ── Helpers ───────────────────────────────────────────────────────────────
    bool is_optimal()  const { return status == SolverStatus::OPTIMAL; }
    bool is_feasible() const {
        return status == SolverStatus::OPTIMAL ||
               status == SolverStatus::TIME_LIMIT ||
               status == SolverStatus::ITERATION_LIMIT;
    }
};

} // namespace suplex
