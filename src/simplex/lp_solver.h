#pragma once
#include "../core/types.h"
#include "../core/problem.h"
#include "../core/solution.h"
#include <vector>
#include <memory>

namespace suplex {

// ── LP Solver interface ───────────────────────────────────────────────────────
///
/// Abstract base used by Person 3's MILPSolver to solve LP relaxations.
/// Person 1 implements PrimalSimplexSolver and DualSimplexSolver.
/// Person 2 implements IPMSolver.
///
class LPSolver {
public:
    virtual ~LPSolver() = default;

    /// Solve from scratch.
    virtual SolverStatus solve(const Problem& problem, Solution& solution) = 0;

    /// Warm-start solve. `basis_indices[i]` = column of i-th basic variable.
    /// Called thousands of times during B&B — MUST be fast.
    virtual SolverStatus solve_from_basis(
        const Problem&          problem,
        const std::vector<Index>& basis_indices,
        Solution&               solution) = 0;

    virtual void set_iteration_limit(Index limit) = 0;
    virtual void set_time_limit(Real seconds)      = 0;
    virtual void set_log_level(LogLevel level)     = 0;

    /// Current basis after solve — used by B&B to warm-start child nodes.
    virtual const std::vector<Index>& get_basis() const = 0;

    /// Farkas infeasibility certificate (used by conflict analysis).
    virtual const std::vector<Real>& get_farkas_ray() const = 0;
};

} // namespace suplex
