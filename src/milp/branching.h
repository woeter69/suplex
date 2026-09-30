#pragma once
#include "../core/types.h"
#include "../core/problem.h"
#include "../core/solution.h"
#include "../simplex/lp_solver.h"
#include "node.h"
#include <vector>
#include <cmath>

namespace suplex {

// ── Branching variable selection ──────────────────────────────────────────────

/// Per-variable pseudo-cost statistics, accumulated across all B&B nodes.
struct PseudoCosts {
    std::vector<Real>  down_cost;   ///< Average obj change per unit when branching down
    std::vector<Real>  up_cost;     ///< Average obj change per unit when branching up
    std::vector<Index> down_count;  ///< Number of down observations
    std::vector<Index> up_count;    ///< Number of up observations

    explicit PseudoCosts(Index num_cols);

    void update_down(Index j, Real frac_val, Real delta_obj);
    void update_up  (Index j, Real frac_val, Real delta_obj);

    bool is_reliable(Index j, Index threshold = 8) const;
    Real score(Index j, Real frac_val) const;  ///< Product score
};

// ── Branching engine ──────────────────────────────────────────────────────────
class Branching {
public:
    explicit Branching(Index num_cols);

    /// Select the best variable to branch on.
    /// Returns the column index, or -1 if all integer vars are integral.
    Index select_variable(const Problem&  p,
                          const Solution& lp_sol,
                          LPSolver*       lp_solver,   ///< For strong branching
                          Real            current_bound,
                          int             node_depth);

    /// Called after a B&B node is solved to record pseudo-cost observations.
    void record_observation(Index j, BranchDirection dir,
                            Real frac_val,
                            Real parent_bound, Real child_bound);

    void set_reliability_threshold(Index t) { reliability_threshold_ = t; }
    void set_max_strong_branch_candidates(Index k) { strong_branch_k_ = k; }

private:
    PseudoCosts pseudo_costs_;
    Index reliability_threshold_  = 8;
    Index strong_branch_k_        = 20;   ///< Top-K candidates for strong branching
    Index strong_branch_max_depth_= 5;    ///< Only strong-branch at shallow depths

    // ── Scoring strategies ────────────────────────────────────────────────
    Index most_fractional      (const Problem& p, const Solution& lp_sol) const;
    Index pseudo_cost_branch   (const Problem& p, const Solution& lp_sol) const;
    Index strong_branch        (const Problem& p, const Solution& lp_sol,
                                LPSolver* lp_solver, Real current_bound);
    Index reliability_branch   (const Problem& p, const Solution& lp_sol,
                                LPSolver* lp_solver, Real current_bound,
                                int node_depth);

    Real  fractional_part(Real val) const { return val - std::floor(val); }
};

} // namespace suplex
