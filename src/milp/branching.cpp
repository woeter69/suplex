#include "branching.h"
#include "../simplex/lp_solver.h"
#include <algorithm>
#include <cmath>
#include <cassert>
#include <numeric>

namespace suplex {

// ── PseudoCosts ───────────────────────────────────────────────────────────────

PseudoCosts::PseudoCosts(Index num_cols)
    : down_cost(num_cols, 1.0),
      up_cost  (num_cols, 1.0),
      down_count(num_cols, 0),
      up_count  (num_cols, 0) {}

void PseudoCosts::update_down(Index j, Real frac_val, Real delta_obj) {
    Real f = frac_val - std::floor(frac_val);
    if (f < 1e-10) return;
    Real obs = delta_obj / f;
    // Running average
    down_cost[j] = (down_cost[j] * down_count[j] + obs) / (down_count[j] + 1);
    down_count[j]++;
}

void PseudoCosts::update_up(Index j, Real frac_val, Real delta_obj) {
    Real f = std::ceil(frac_val) - frac_val;
    if (f < 1e-10) return;
    Real obs = delta_obj / f;
    up_cost[j] = (up_cost[j] * up_count[j] + obs) / (up_count[j] + 1);
    up_count[j]++;
}

bool PseudoCosts::is_reliable(Index j, Index threshold) const {
    return down_count[j] >= threshold && up_count[j] >= threshold;
}

/// Product scoring (SCIP default, μ = 1/6)
Real PseudoCosts::score(Index j, Real frac_val) const {
    Real f      = frac_val - std::floor(frac_val);
    Real d_est  = down_cost[j] * f;
    Real u_est  = up_cost[j]   * (1.0 - f);
    constexpr Real mu = 1.0 / 6.0;
    return (1.0 - mu) * std::min(d_est, u_est) + mu * std::max(d_est, u_est);
}

// ── Branching ─────────────────────────────────────────────────────────────────

Branching::Branching(Index num_cols)
    : pseudo_costs_(num_cols) {}

void Branching::record_observation(Index j, BranchDirection dir,
                                   Real frac_val,
                                   Real parent_bound, Real child_bound) {
    Real delta = child_bound - parent_bound;
    if (delta < 0.0) delta = 0.0;   // bound cannot decrease
    if (dir == BranchDirection::DOWN)
        pseudo_costs_.update_down(j, frac_val, delta);
    else
        pseudo_costs_.update_up  (j, frac_val, delta);
}

Index Branching::select_variable(const Problem&  p,
                                  const Solution& lp_sol,
                                  LPSolver*       lp_solver,
                                  Real            current_bound,
                                  int             node_depth) {
    return reliability_branch(p, lp_sol, lp_solver, current_bound, node_depth);
}

// ── Most fractional ────────────────────────────────────────────────────────────

Index Branching::most_fractional(const Problem& p, const Solution& lp_sol) const {
    const auto& vt = p.var_types();
    Index best_j = -1;
    Real  best_f = -1.0;

    for (Index j = 0; j < p.num_cols(); ++j) {
        if (vt[j] == VarType::CONTINUOUS) continue;
        Real val = lp_sol.primal_values[j];
        Real f   = val - std::floor(val);   // fractional part in [0,1)
        Real dist = std::min(f, 1.0 - f);   // distance to nearest integer
        if (dist > EPS_INTEGER && dist > best_f) {
            best_f = dist;
            best_j = j;
        }
    }
    return best_j;
}

// ── Pseudo-cost branching ──────────────────────────────────────────────────────

Index Branching::pseudo_cost_branch(const Problem& p, const Solution& lp_sol) const {
    const auto& vt = p.var_types();
    Index best_j    = -1;
    Real  best_score = -1.0;

    for (Index j = 0; j < p.num_cols(); ++j) {
        if (vt[j] == VarType::CONTINUOUS) continue;
        Real val  = lp_sol.primal_values[j];
        Real frac = val - std::floor(val);
        if (frac < EPS_INTEGER || frac > 1.0 - EPS_INTEGER) continue;  // already integer

        Real sc = pseudo_costs_.score(j, val);
        if (sc > best_score) {
            best_score = sc;
            best_j     = j;
        }
    }
    return best_j;
}

// ── Strong branching ───────────────────────────────────────────────────────────

Index Branching::strong_branch(const Problem& p, const Solution& lp_sol,
                                LPSolver* lp_solver, Real current_bound) {
    if (!lp_solver) return most_fractional(p, lp_sol);

    const auto& vt = p.var_types();
    const auto& cl = p.col_lower();
    const auto& cu = p.col_upper();

    // Collect fractional variables, sort by fractionality (closest to 0.5 first)
    std::vector<std::pair<Real, Index>> candidates;
    for (Index j = 0; j < p.num_cols(); ++j) {
        if (vt[j] == VarType::CONTINUOUS) continue;
        Real val  = lp_sol.primal_values[j];
        Real frac = val - std::floor(val);
        if (frac < EPS_INTEGER || frac > 1.0 - EPS_INTEGER) continue;
        candidates.emplace_back(std::min(frac, 1.0 - frac), j);
    }
    if (candidates.empty()) return -1;

    // Sort descending by fractionality (most fractional = closest to 0.5)
    std::sort(candidates.begin(), candidates.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    Index K = std::min(strong_branch_k_, static_cast<Index>(candidates.size()));

    Index best_j    = candidates[0].second;
    Real  best_score = -1.0;

    for (Index ci = 0; ci < K; ++ci) {
        Index j   = candidates[ci].second;
        Real  val = lp_sol.primal_values[j];

        // Try branching DOWN: x_j <= floor(val)
        Problem p_down = p;
        {
            auto new_cu = cu;
            new_cu[j]   = std::floor(val);
            p_down.set_col_bounds(std::vector<Real>(cl), std::move(new_cu));
        }
        Solution sol_down;
        lp_solver->solve(p_down, sol_down);
        Real down_bound = sol_down.is_optimal()
                         ? sol_down.objective_value : 1e30;

        // Try branching UP: x_j >= ceil(val)
        Problem p_up = p;
        {
            auto new_cl = cl;
            new_cl[j]   = std::ceil(val);
            p_up.set_col_bounds(std::move(new_cl), std::vector<Real>(cu));
        }
        Solution sol_up;
        lp_solver->solve(p_up, sol_up);
        Real up_bound = sol_up.is_optimal()
                       ? sol_up.objective_value : 1e30;

        // Record pseudo-cost observations from strong branching
        if (sol_down.is_optimal())
            pseudo_costs_.update_down(j, val, down_bound - current_bound);
        if (sol_up.is_optimal())
            pseudo_costs_.update_up  (j, val, up_bound   - current_bound);

        Real sc = pseudo_costs_.score(j, val);
        if (sc > best_score) {
            best_score = sc;
            best_j     = j;
        }
    }
    return best_j;
}

// ── Reliability branching (default) ───────────────────────────────────────────

Index Branching::reliability_branch(const Problem& p, const Solution& lp_sol,
                                    LPSolver* lp_solver, Real current_bound,
                                    int node_depth) {
    const auto& vt = p.var_types();

    // Classify fractional variables
    std::vector<Index> reliable, unreliable;
    for (Index j = 0; j < p.num_cols(); ++j) {
        if (vt[j] == VarType::CONTINUOUS) continue;
        Real val  = lp_sol.primal_values[j];
        Real frac = val - std::floor(val);
        if (frac < EPS_INTEGER || frac > 1.0 - EPS_INTEGER) continue;

        if (pseudo_costs_.is_reliable(j, reliability_threshold_))
            reliable.push_back(j);
        else
            unreliable.push_back(j);
    }

    if (reliable.empty() && unreliable.empty()) return -1;  // all integer

    // Always strong-branch on unreliable candidates at shallow depths
    bool do_strong = lp_solver &&
                     !unreliable.empty() &&
                     node_depth <= strong_branch_max_depth_;

    if (do_strong) {
        // Sort unreliable by fractionality (closest to 0.5)
        std::sort(unreliable.begin(), unreliable.end(), [&](Index a, Index b) {
            Real fa = lp_sol.primal_values[a] - std::floor(lp_sol.primal_values[a]);
            Real fb = lp_sol.primal_values[b] - std::floor(lp_sol.primal_values[b]);
            return std::min(fa, 1-fa) > std::min(fb, 1-fb);
        });

        Index K = std::min(strong_branch_k_, static_cast<Index>(unreliable.size()));
        unreliable.resize(K);

        const auto& cl = p.col_lower();
        const auto& cu = p.col_upper();

        Index best_j    = unreliable[0];
        Real  best_score = -1.0;

        for (Index j : unreliable) {
            Real val = lp_sol.primal_values[j];

            // Down
            Problem p_down = p;
            auto new_cu = cu; new_cu[j] = std::floor(val);
            p_down.set_col_bounds(std::vector<Real>(cl), std::move(new_cu));
            Solution sd; lp_solver->solve(p_down, sd);
            Real db = sd.is_optimal() ? sd.objective_value : 1e30;
            if (sd.is_optimal()) pseudo_costs_.update_down(j, val, db - current_bound);

            // Up
            Problem p_up = p;
            auto new_cl = cl; new_cl[j] = std::ceil(val);
            p_up.set_col_bounds(std::move(new_cl), std::vector<Real>(cu));
            Solution su; lp_solver->solve(p_up, su);
            Real ub = su.is_optimal() ? su.objective_value : 1e30;
            if (su.is_optimal()) pseudo_costs_.update_up(j, val, ub - current_bound);

            Real sc = pseudo_costs_.score(j, val);
            if (sc > best_score) { best_score = sc; best_j = j; }
        }
        return best_j;
    }

    // Use pseudo-cost branching for reliable variables (or fallback)
    if (!reliable.empty()) return pseudo_cost_branch(p, lp_sol);
    return most_fractional(p, lp_sol);
}

} // namespace suplex
