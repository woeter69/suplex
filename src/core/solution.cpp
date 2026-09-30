#include "solution.h"

namespace suplex {

bool Solution::is_feasible() const {
    return status == SolverStatus::OPTIMAL ||
           status == SolverStatus::ITERATION_LIMIT ||
           status == SolverStatus::TIME_LIMIT;
}

void Solution::reset() {
    status = SolverStatus::NOT_STARTED;
    objective_value = 0.0;
    primal_values.clear();
    dual_values.clear();
    reduced_costs.clear();
    row_activities.clear();
    best_bound = 0.0;
    mip_gap = 0.0;
    num_nodes = 0;
    num_iterations = 0;
    solve_time_seconds = 0.0;
}

} // namespace suplex
