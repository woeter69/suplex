#pragma once
#include "../core/types.h"
#include "../core/problem.h"
#include "../core/solution.h"
#include "presolve.h"
#include "branch_bound.h"   // MILPCallback, MILPSolver, TreeStats
#include "../simplex/lp_solver.h"
#include <memory>
#include <vector>

namespace suplex {

// ── MILP Solver public API ────────────────────────────────────────────────────
///
/// Top-level solver for mixed-integer programs.
/// Wires together: Presolve → Root LP → Cuts → Heuristics → B&B → Postsolve.
///
class MILPSolver {
public:
    MILPSolver();

    /// Set the LP solver used to solve relaxations (Person 1/2 provide this).
    void set_lp_solver(std::shared_ptr<LPSolver> lp_solver);

    /// Solve. Problem must already be set (or pass directly).
    SolverStatus solve(const Problem& problem, Solution& solution);

    // ── Configuration ─────────────────────────────────────────────────────
    void set_gap_tolerance  (Real gap)          { gap_tol_     = gap; }
    void set_node_limit     (int64_t limit)     { node_limit_  = limit; }
    void set_time_limit     (Real seconds)      { time_limit_  = seconds; }
    void set_log_level      (LogLevel level)    { log_level_   = level; }
    void set_presolve       (bool enabled)      { use_presolve_= enabled; }
    void set_callback       (std::shared_ptr<MILPCallback> cb);

    /// Access solve statistics.
    const TreeStats& tree_stats() const;

private:
    std::shared_ptr<LPSolver>     lp_solver_;
    std::shared_ptr<MILPCallback> callback_;

    Real     gap_tol_     = 1e-4;
    int64_t  node_limit_  = 1000000000LL;
    Real     time_limit_  = 3600.0;
    LogLevel log_level_   = LogLevel::INFO;
    bool     use_presolve_= true;
};

} // namespace suplex
