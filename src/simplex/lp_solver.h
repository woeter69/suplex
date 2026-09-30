#pragma once

#include "src/core/types.h"
#include "src/core/problem.h"
#include "src/core/solution.h"
#include <vector>

namespace suplex {

/**
 * @brief Abstract interface for Linear Programming solvers.
 * Implemented by PrimalSimplexSolver, DualSimplexSolver, and IPMSolver.
 */
class LPSolver {
public:
    virtual ~LPSolver() = default;

    /**
     * @brief Solves the given linear programming problem from scratch.
     */
    virtual SolverStatus solve(const Problem& problem, Solution& solution) = 0;

    /**
     * @brief Warm-starts the solve from a provided basis.
     * Essential for MILP branch-and-bound iterations.
     */
    virtual SolverStatus solve_from_basis(
        const Problem& problem,
        const std::vector<Index>& basis_indices,
        Solution& solution
    ) = 0;

    // Configuration
    virtual void set_iteration_limit(Index limit) = 0;
    virtual void set_time_limit(Real seconds) = 0;
    virtual void set_log_level(LogLevel level) = 0;

    // Basis inspection
    virtual const std::vector<Index>& get_basis() const = 0;
};

} // namespace suplex
