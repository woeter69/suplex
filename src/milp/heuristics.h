#pragma once
#include "../core/types.h"
#include "../core/problem.h"
#include "../core/solution.h"
#include "../simplex/lp_solver.h"
#include <vector>
#include <optional>

namespace suplex {

// ── Primal heuristics ─────────────────────────────────────────────────────────
///
/// All heuristics accept an LP relaxation solution and try to find an
/// integer-feasible solution.  Return std::nullopt if they fail.
///
class Heuristics {
public:
    explicit Heuristics(LPSolver* lp_solver);

    // ── Simple rounding ───────────────────────────────────────────────────
    /// Round each fractional integer variable to nearest integer.
    /// O(n) — fast but low success rate.
    std::optional<Solution> simple_rounding(const Problem&  p,
                                            const Solution& lp_sol);

    // ── Fractional diving ─────────────────────────────────────────────────
    /// Repeatedly fix the most fractional variable and re-solve via dual simplex.
    std::optional<Solution> fractional_diving(const Problem&  p,
                                              const Solution& lp_sol,
                                              int max_dives = 100);

    // ── Coefficient diving ────────────────────────────────────────────────
    /// Fix variables that have the smallest objective impact.
    std::optional<Solution> coefficient_diving(const Problem&  p,
                                               const Solution& lp_sol,
                                               int max_dives = 100);

    // ── Feasibility pump ──────────────────────────────────────────────────
    /// Alternate between rounding and projecting back to LP feasibility.
    std::optional<Solution> feasibility_pump(const Problem&  p,
                                             const Solution& lp_sol,
                                             int max_iter = 100);

    // ── RINS ─────────────────────────────────────────────────────────────
    /// Fix agreement variables, solve sub-MIP with node limit.
    std::optional<Solution> rins(const Problem&  p,
                                 const Solution& lp_sol,
                                 const Solution& incumbent,
                                 int64_t max_nodes = 500);

    void set_log_level(LogLevel l) { log_level_ = l; }

private:
    LPSolver* lp_solver_;
    LogLevel  log_level_ = LogLevel::INFO;

    // ── Helpers ───────────────────────────────────────────────────────────
    bool is_feasible         (const Problem& p, const Solution& sol) const;
    bool is_integer_feasible (const Problem& p, const Solution& sol) const;
    void round_integers      (const Problem& p, Solution& sol) const;

    /// Solve re-optimisation LP with modified bounds; returns false if infeasible.
    bool resolve_lp(const Problem& p,
                    const std::vector<Real>& temp_col_lower,
                    const std::vector<Real>& temp_col_upper,
                    Solution& sol);
};

} // namespace suplex
