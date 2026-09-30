#include "heuristics.h"
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>

namespace suplex {

Heuristics::Heuristics(LPSolver* lp_solver) : lp_solver_(lp_solver) {}

// ── Helper: check feasibility ─────────────────────────────────────────────────

bool Heuristics::is_feasible(const Problem& p, const Solution& sol) const {
    const auto& rl = p.row_lower();
    const auto& ru = p.row_upper();
    const auto& cl = p.col_lower();
    const auto& cu = p.col_upper();
    const auto& x  = sol.primal_values;

    for (Index j = 0; j < p.num_cols(); ++j)
        if (x[j] < cl[j] - EPS_FEASIBILITY || x[j] > cu[j] + EPS_FEASIBILITY)
            return false;

    // Check row activities (need Ax)
    std::vector<Real> ax(p.num_rows(), 0.0);
    p.constraint_matrix().multiply_vec(x.data(), ax.data());
    for (Index i = 0; i < p.num_rows(); ++i)
        if (ax[i] < rl[i] - EPS_FEASIBILITY || ax[i] > ru[i] + EPS_FEASIBILITY)
            return false;

    return true;
}

bool Heuristics::is_integer_feasible(const Problem& p, const Solution& sol) const {
    const auto& vt = p.var_types();
    for (Index j = 0; j < p.num_cols(); ++j) {
        if (vt[j] == VarType::CONTINUOUS) continue;
        Real f = sol.primal_values[j] - std::floor(sol.primal_values[j]);
        if (f > EPS_INTEGER && f < 1.0 - EPS_INTEGER) return false;
    }
    return true;
}

void Heuristics::round_integers(const Problem& p, Solution& sol) const {
    const auto& vt = p.var_types();
    const auto& cl = p.col_lower();
    const auto& cu = p.col_upper();
    for (Index j = 0; j < p.num_cols(); ++j) {
        if (vt[j] == VarType::CONTINUOUS) continue;
        Real val = sol.primal_values[j];
        Real rounded = std::round(val);
        // Clamp to bounds
        sol.primal_values[j] = std::max(cl[j], std::min(cu[j], rounded));
    }
}

bool Heuristics::resolve_lp(const Problem& p,
                              const std::vector<Real>& temp_cl,
                              const std::vector<Real>& temp_cu,
                              Solution& sol) {
    Problem p_mod = p;
    p_mod.set_col_bounds(std::vector<Real>(temp_cl), std::vector<Real>(temp_cu));

    SolverStatus status = lp_solver_->solve(p_mod, sol);
    return status == SolverStatus::OPTIMAL;
}

// ── 1. Simple rounding ────────────────────────────────────────────────────────

std::optional<Solution> Heuristics::simple_rounding(const Problem&  p,
                                                      const Solution& lp_sol) {
    Solution candidate = lp_sol;
    round_integers(p, candidate);

    if (!is_feasible(p, candidate)) return std::nullopt;

    // Compute objective value
    const auto& c = p.objective();
    candidate.objective_value = 0.0;
    for (Index j = 0; j < p.num_cols(); ++j)
        candidate.objective_value += c[j] * candidate.primal_values[j];
    candidate.status = SolverStatus::OPTIMAL;

    if (log_level_ <= LogLevel::DEBUG)
        std::cout << "[Heuristic] Simple rounding: obj=" << candidate.objective_value << "\n";

    return candidate;
}

// ── 2. Fractional diving ──────────────────────────────────────────────────────

std::optional<Solution> Heuristics::fractional_diving(const Problem&  p,
                                                        const Solution& lp_sol,
                                                        int max_dives) {
    if (!lp_solver_) return std::nullopt;

    std::vector<Real> cur_cl = p.col_lower();
    std::vector<Real> cur_cu = p.col_upper();
    Solution cur_sol = lp_sol;

    const auto& vt = p.var_types();

    for (int dive = 0; dive < max_dives; ++dive) {
        // Find most fractional integer variable
        Index best_j  = -1;
        Real  best_f  = -1.0;

        for (Index j = 0; j < p.num_cols(); ++j) {
            if (vt[j] == VarType::CONTINUOUS) continue;
            Real val  = cur_sol.primal_values[j];
            Real frac = val - std::floor(val);
            Real dist = std::min(frac, 1.0 - frac);
            if (dist > EPS_INTEGER && dist > best_f) {
                best_f = dist;
                best_j = j;
            }
        }

        if (best_j == -1) {
            // All integer — check feasibility
            if (is_feasible(p, cur_sol)) {
                const auto& c = p.objective();
                cur_sol.objective_value = 0.0;
                for (Index j = 0; j < p.num_cols(); ++j)
                    cur_sol.objective_value += c[j] * cur_sol.primal_values[j];
                cur_sol.status = SolverStatus::OPTIMAL;
                if (log_level_ <= LogLevel::DEBUG)
                    std::cout << "[Heuristic] Fractional diving: obj=" << cur_sol.objective_value << "\n";
                return cur_sol;
            }
            return std::nullopt;
        }

        // Fix to nearest integer
        Real val     = cur_sol.primal_values[best_j];
        Real rounded = (val - std::floor(val) <= 0.5)
                       ? std::floor(val) : std::ceil(val);
        rounded = std::max(cur_cl[best_j], std::min(cur_cu[best_j], rounded));

        // Tighten bound (both lb and ub = rounded, effectively fixing it)
        Real saved_lb = cur_cl[best_j];
        Real saved_ub = cur_cu[best_j];
        cur_cl[best_j] = rounded;
        cur_cu[best_j] = rounded;

        bool ok = resolve_lp(p, cur_cl, cur_cu, cur_sol);
        if (!ok) {
            // LP infeasible — dive failed
            cur_cl[best_j] = saved_lb;
            cur_cu[best_j] = saved_ub;
            return std::nullopt;
        }
    }
    return std::nullopt;
}

// ── 3. Coefficient diving ─────────────────────────────────────────────────────

std::optional<Solution> Heuristics::coefficient_diving(const Problem&  p,
                                                         const Solution& lp_sol,
                                                         int max_dives) {
    if (!lp_solver_) return std::nullopt;

    std::vector<Real> cur_cl = p.col_lower();
    std::vector<Real> cur_cu = p.col_upper();
    Solution cur_sol = lp_sol;
    const auto& vt = p.var_types();
    const auto& c  = p.objective();

    for (int dive = 0; dive < max_dives; ++dive) {
        // Find fractional var with smallest |c_j| (least objective impact)
        Index best_j     = -1;
        Real  best_score = 1e30;

        for (Index j = 0; j < p.num_cols(); ++j) {
            if (vt[j] == VarType::CONTINUOUS) continue;
            Real val  = cur_sol.primal_values[j];
            Real frac = val - std::floor(val);
            if (frac < EPS_INTEGER || frac > 1.0 - EPS_INTEGER) continue;
            Real score = std::abs(c[j]);
            if (score < best_score) { best_score = score; best_j = j; }
        }

        if (best_j == -1) {
            if (is_feasible(p, cur_sol)) {
                cur_sol.objective_value = 0.0;
                for (Index j = 0; j < p.num_cols(); ++j)
                    cur_sol.objective_value += c[j] * cur_sol.primal_values[j];
                cur_sol.status = SolverStatus::OPTIMAL;
                return cur_sol;
            }
            return std::nullopt;
        }

        Real val     = cur_sol.primal_values[best_j];
        Real rounded = (c[best_j] >= 0) ? std::floor(val) : std::ceil(val);
        rounded = std::max(cur_cl[best_j], std::min(cur_cu[best_j], rounded));

        cur_cl[best_j] = rounded;
        cur_cu[best_j] = rounded;

        if (!resolve_lp(p, cur_cl, cur_cu, cur_sol)) return std::nullopt;
    }
    return std::nullopt;
}

// ── 4. Feasibility pump ───────────────────────────────────────────────────────

std::optional<Solution> Heuristics::feasibility_pump(const Problem&  p,
                                                       const Solution& lp_sol,
                                                       int max_iter) {
    if (!lp_solver_) return std::nullopt;

    const auto& vt = p.var_types();
    const auto& cl = p.col_lower();
    const auto& cu = p.col_upper();
    const auto& c  = p.objective();

    Solution x_star = lp_sol;  // current LP solution
    std::vector<Real> x_tilde(p.num_cols());   // rounded target
    std::vector<Real> prev_tilde(p.num_cols(), -1e30);

    Real alpha = 0.85;   // weight on distance vs. original objective

    for (int iter = 0; iter < max_iter; ++iter) {
        // Round x_star to x_tilde
        for (Index j = 0; j < p.num_cols(); ++j) {
            if (vt[j] == VarType::CONTINUOUS)
                x_tilde[j] = x_star.primal_values[j];
            else {
                Real rounded = std::round(x_star.primal_values[j]);
                x_tilde[j] = std::max(cl[j], std::min(cu[j], rounded));
            }
        }

        // Check if x_tilde is feasible
        Solution candidate;
        candidate.primal_values = x_tilde;
        candidate.status = SolverStatus::OPTIMAL;
        candidate.objective_value = 0.0;
        for (Index j = 0; j < p.num_cols(); ++j)
            candidate.objective_value += c[j] * x_tilde[j];

        if (is_feasible(p, candidate) && is_integer_feasible(p, candidate)) {
            if (log_level_ <= LogLevel::DEBUG)
                std::cout << "[Heuristic] Feasibility pump iter=" << iter
                          << " obj=" << candidate.objective_value << "\n";
            return candidate;
        }

        // Detect cycling: if x_tilde == prev_tilde, perturb
        bool cycling = true;
        for (Index j = 0; j < p.num_cols() && cycling; ++j)
            if (std::abs(x_tilde[j] - prev_tilde[j]) > 0.5) cycling = false;

        if (cycling) {
            // Randomly flip some binary roundings
            for (Index j = 0; j < p.num_cols(); ++j) {
                if (vt[j] != VarType::BINARY && vt[j] != VarType::INTEGER) continue;
                // ~10% perturbation probability
                if ((std::rand() % 10) == 0) {
                    x_tilde[j] = (x_tilde[j] == std::floor(x_tilde[j]))
                                 ? std::ceil(x_tilde[j]) : std::floor(x_tilde[j]);
                    x_tilde[j] = std::max(cl[j], std::min(cu[j], x_tilde[j]));
                }
            }
        }
        prev_tilde = x_tilde;

        // Solve LP: min alpha * Σ|x_j - x_tilde_j| + (1-alpha) * c^T x
        // Approximated by: modify objective c_j' = (1-alpha)*c_j ± alpha * sign
        // (full L1 implementation uses auxiliary variables — simplified here)
        std::vector<Real> c_pump(p.num_cols());
        for (Index j = 0; j < p.num_cols(); ++j) {
            Real diff = x_star.primal_values[j] - x_tilde[j];
            Real sign = (diff > 0) ? 1.0 : -1.0;
            if (vt[j] == VarType::CONTINUOUS)
                c_pump[j] = (1.0 - alpha) * c[j];
            else
                c_pump[j] = (1.0 - alpha) * c[j] + alpha * sign;
        }

        // Build modified problem with pump objective
        Problem p_pump = p;
        p_pump.set_objective(std::vector<Real>(c_pump));

        SolverStatus st = lp_solver_->solve(p_pump, x_star);
        if (st != SolverStatus::OPTIMAL) break;

        alpha *= 0.95;  // gradually reduce weight on distance
    }

    return std::nullopt;
}

// ── 5. RINS ───────────────────────────────────────────────────────────────────

std::optional<Solution> Heuristics::rins(const Problem&  p,
                                          const Solution& lp_sol,
                                          const Solution& incumbent,
                                          int64_t         max_nodes) {
    if (!lp_solver_) return std::nullopt;

    const auto& vt = p.var_types();
    const auto& cl = p.col_lower();
    const auto& cu = p.col_upper();

    // Fix variables where LP relaxation and incumbent agree
    std::vector<Real> fixed_cl = cl;
    std::vector<Real> fixed_cu = cu;
    Index num_fixed = 0;

    for (Index j = 0; j < p.num_cols(); ++j) {
        if (vt[j] == VarType::CONTINUOUS) continue;
        Real lp_val  = lp_sol.primal_values[j];
        Real inc_val = incumbent.primal_values[j];
        Real lp_int  = std::round(lp_val);
        Real inc_int = std::round(inc_val);

        // Both agree on the same integer value
        if (std::abs(lp_int - inc_int) < 0.5 &&
            lp_int >= cl[j] - EPS_FEASIBILITY &&
            lp_int <= cu[j] + EPS_FEASIBILITY) {
            fixed_cl[j] = lp_int;
            fixed_cu[j] = lp_int;
            num_fixed++;
        }
    }

    if (log_level_ <= LogLevel::DEBUG)
        std::cout << "[RINS] Fixing " << num_fixed << "/" << p.num_cols() << " variables\n";

    // Solve restricted LP (not a full MIP solve — would need recursive call)
    // For now, solve the LP relaxation of the fixed problem
    Problem p_fixed = p;
    p_fixed.set_col_bounds(std::move(fixed_cl), std::move(fixed_cu));
    Solution sol;
    SolverStatus st = lp_solver_->solve(p_fixed, sol);
    if (st != SolverStatus::OPTIMAL) return std::nullopt;

    // Round the unfixed variables
    round_integers(p_fixed, sol);
    if (!is_feasible(p, sol)) return std::nullopt;
    if (!is_integer_feasible(p, sol)) return std::nullopt;

    const auto& c = p.objective();
    sol.objective_value = 0.0;
    for (Index j = 0; j < p.num_cols(); ++j)
        sol.objective_value += c[j] * sol.primal_values[j];
    sol.status = SolverStatus::OPTIMAL;
    return sol;
}

} // namespace suplex
