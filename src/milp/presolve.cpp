#include "presolve.h"
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>
#include <numeric>

namespace suplex {

// ── Helpers ──────────────────────────────────────────────────────────────────

static inline bool nearly_equal(Real a, Real b, Real tol = EPS_FEASIBILITY) {
    return std::abs(a - b) <= tol;
}
static inline bool is_inf(Real v) { return v >= INF || v <= -INF; }

// ── Public API ────────────────────────────────────────────────────────────────

PresolveResult Presolve::apply(const Problem& original) {
    PresolveResult result;
    result.obj_offset = 0.0;

    // Work on a mutable copy
    Problem p = original;

    bool infeasible = false;

    if (log_level_ <= LogLevel::INFO) {
        std::cout << "[Presolve] Starting: "
                  << p.num_rows() << " rows, "
                  << p.num_cols() << " cols\n";
    }

    for (int round = 0; round < max_rounds_; ++round) {
        bool changed = false;

        changed |= reduce_empty_rows        (p, result.stack, infeasible);
        if (infeasible) { result.is_infeasible = true; return result; }

        changed |= reduce_fixed_variables   (p, result.stack);
        changed |= reduce_singleton_rows    (p, result.stack, infeasible);
        if (infeasible) { result.is_infeasible = true; return result; }

        changed |= reduce_singleton_cols    (p, result.stack);
        changed |= reduce_forcing_rows      (p, result.stack, infeasible);
        if (infeasible) { result.is_infeasible = true; return result; }

        changed |= tighten_implied_bounds   (p, result.stack, infeasible);
        if (infeasible) { result.is_infeasible = true; return result; }

        changed |= reduce_duplicate_rows    (p, result.stack);
        changed |= reduce_duplicate_cols    (p, result.stack);

        if (p.is_mip()) {
            changed |= tighten_coefficients (p, result.stack);
            if (enable_probing_) {
                changed |= run_probing      (p, result.stack, infeasible);
                if (infeasible) { result.is_infeasible = true; return result; }
            }
        }

        if (!changed) break;
    }

    if (log_level_ <= LogLevel::INFO) {
        std::cout << "[Presolve] Done: "
                  << p.num_rows() << " rows, "
                  << p.num_cols() << " cols, "
                  << "offset=" << result.obj_offset << "\n";
    }

    result.reduced_problem = std::move(p);
    return result;
}

// ── Postsolve ─────────────────────────────────────────────────────────────────

Solution Presolve::postsolve(const Solution& reduced_sol,
                              const PresolveResult& result) {
    // Start with the reduced-space solution and expand it
    // Work on a mutable copy
    Solution sol = reduced_sol;

    // We need to map reduced indices → original indices.
    // The stack was built in apply() order; pop in reverse.
    PresolveStack& stack = const_cast<PresolveStack&>(result.stack);

    // Expand primal, dual, rc to original dimensions
    // (In a full implementation we'd track the mapping; for now we
    //  process each record and adjust in-place.)

    while (!stack.empty()) {
        PresolveRecord rec = stack.pop();

        switch (rec.rule) {
        case PresolveRuleType::FIXED_VARIABLE: {
            // x_j was fixed and removed. Restore it.
            Index j = rec.col_idx;
            if (j < static_cast<Index>(sol.primal_values.size()))
                sol.primal_values[j] = rec.fix_value;
            // reduced cost: c_j - y^T a_j (approximately)
            break;
        }
        case PresolveRuleType::SINGLETON_ROW: {
            // Row was removed; variable bounds were tightened.
            // Row dual y[i] depends on which bound is active.
            Index i = rec.row_idx;
            Index j = rec.col_idx;
            Real xj = (j < static_cast<Index>(sol.primal_values.size()))
                      ? sol.primal_values[j] : rec.saved_lb;

            Real y_i = 0.0;
            // If constraint is tight at row_upper (original): y_i >= 0
            // If constraint is tight at row_lower (original): y_i <= 0
            if (nearly_equal(rec.a_ij * xj, rec.saved_ru))
                y_i = (sol.reduced_costs.size() > (size_t)j)
                      ? sol.reduced_costs[j] / rec.a_ij : 0.0;
            if (i < static_cast<Index>(sol.dual_values.size()))
                sol.dual_values[i] = y_i;
            break;
        }
        case PresolveRuleType::SINGLETON_COL: {
            // Restore x_j = optimal value
            Index j = rec.col_idx;
            if (j < static_cast<Index>(sol.primal_values.size()))
                sol.primal_values[j] = rec.fix_value;
            // Row dual: y_i = c_j / a_ij if x_j strictly between bounds
            Index i = rec.row_idx;
            if (i < static_cast<Index>(sol.dual_values.size()) &&
                std::abs(rec.a_ij) > EPS_ZERO)
                sol.dual_values[i] = rec.c_j / rec.a_ij;
            break;
        }
        case PresolveRuleType::FORCING_ROW: {
            // All variables fixed at bounds; restore them
            for (Index k = 0;
                 k < static_cast<Index>(rec.affected_cols.size()); ++k) {
                Index j = rec.affected_cols[k];
                if (j < static_cast<Index>(sol.primal_values.size()))
                    sol.primal_values[j] = rec.saved_col_lbs[k];
            }
            break;
        }
        case PresolveRuleType::PROBING:
        case PresolveRuleType::IMPLIED_BOUND:
        case PresolveRuleType::COEFFICIENT_TIGHTENING: {
            // Bounds were tightened; no primal value to restore,
            // just re-check feasibility implicitly.
            break;
        }
        default:
            break;
        }
    }

    // Add back objective constant from substitutions
    sol.objective_value += result.obj_offset;
    return sol;
}

// ── Reduction: Empty rows ──────────────────────────────────────────────────────

bool Presolve::reduce_empty_rows(Problem& p, PresolveStack& s, bool& infeasible) {
    // Build CSR view for row iteration
    const auto& A   = p.constraint_matrix();
    const auto& rl  = p.row_lower();
    const auto& ru  = p.row_upper();

    // Count nnz per row
    std::vector<Index> row_nnz(p.num_rows(), 0);
    for (Index j = 0; j < p.num_cols(); ++j)
        for (Index k = A.col_begin(j); k < A.col_end(j); ++k)
            row_nnz[A.row_index[k]]++;

    bool changed = false;
    std::vector<bool> row_active(p.num_rows(), true);

    for (Index i = 0; i < p.num_rows(); ++i) {
        if (row_nnz[i] > 0) continue;
        // Empty row: 0 <= 0 but constraint is rl[i] <= 0 <= ru[i]
        if (rl[i] > EPS_FEASIBILITY || ru[i] < -EPS_FEASIBILITY) {
            // 0 not in [rl, ru] → infeasible
            infeasible = true;
            return true;
        }
        // Row is trivially satisfied; remove it
        row_active[i] = false;
        s.push({PresolveRuleType::EMPTY_ROW, i, -1});
        changed = true;
    }

    if (changed) {
        // Compact problem (remove empty rows)
        std::vector<Index> old_to_new_row, old_to_new_col;
        std::vector<bool>  col_active(p.num_cols(), true);
        p = compact(p, row_active, col_active, old_to_new_row, old_to_new_col);
    }
    return changed;
}

// ── Reduction: Fixed variables ────────────────────────────────────────────────

bool Presolve::reduce_fixed_variables(Problem& p, PresolveStack& s) {
    const auto& cl = p.col_lower();
    const auto& cu = p.col_upper();
    const auto& A  = p.constraint_matrix();

    bool changed = false;
    std::vector<bool> col_active(p.num_cols(), true);
    std::vector<bool> row_active(p.num_rows(), true);

    // We need mutable access to row bounds
    std::vector<Real> new_rl = p.row_lower();
    std::vector<Real> new_ru = p.row_upper();
    std::vector<Real> new_c  = p.objective();
    Real obj_offset = 0.0;

    for (Index j = 0; j < p.num_cols(); ++j) {
        if (!nearly_equal(cl[j], cu[j])) continue;   // not fixed
        Real fix_val = cl[j];

        PresolveRecord rec;
        rec.rule      = PresolveRuleType::FIXED_VARIABLE;
        rec.col_idx   = j;
        rec.fix_value = fix_val;
        rec.saved_lb  = cl[j];
        rec.saved_ub  = cu[j];
        rec.c_j       = new_c[j];
        s.push(rec);

        // Substitute x_j = fix_val into all rows
        for (Index k = A.col_begin(j); k < A.col_end(j); ++k) {
            Index i    = A.row_index[k];
            Real  aij  = A.values[k];
            Real  sub  = aij * fix_val;
            if (!is_inf(new_rl[i])) new_rl[i] -= sub;
            if (!is_inf(new_ru[i])) new_ru[i] -= sub;
        }

        // Update objective constant
        obj_offset += new_c[j] * fix_val;
        new_c[j]    = 0.0;

        col_active[j] = false;
        changed = true;
    }

    if (changed) {
        // Apply updated row bounds and compact
        Problem tmp = p;
        tmp.set_row_bounds(std::move(new_rl), std::move(new_ru));
        tmp.set_objective(std::move(new_c));

        std::vector<Index> m1, m2;
        p = compact(tmp, row_active, col_active, m1, m2);
    }
    return changed;
}

// ── Reduction: Singleton rows ─────────────────────────────────────────────────

bool Presolve::reduce_singleton_rows(Problem& p, PresolveStack& s, bool& infeasible) {
    const auto& A  = p.constraint_matrix();
    const auto& cl = p.col_lower();
    const auto& cu = p.col_upper();
    const auto& rl = p.row_lower();
    const auto& ru = p.row_upper();

    // Count nnz per row
    std::vector<Index> row_nnz(p.num_rows(), 0);
    std::vector<Index> row_col(p.num_rows(), -1);   // column of the singleton
    std::vector<Real>  row_aij(p.num_rows(), 0.0);  // its coefficient

    for (Index j = 0; j < p.num_cols(); ++j)
        for (Index k = A.col_begin(j); k < A.col_end(j); ++k) {
            Index i = A.row_index[k];
            row_nnz[i]++;
            row_col[i] = j;
            row_aij[i] = A.values[k];
        }

    bool changed = false;
    std::vector<bool> row_active(p.num_rows(), true);
    std::vector<Real> new_cl = cl, new_cu = cu;

    for (Index i = 0; i < p.num_rows(); ++i) {
        if (row_nnz[i] != 1) continue;

        Index j   = row_col[i];
        Real  aij = row_aij[i];
        if (std::abs(aij) < EPS_ZERO) continue;

        // Constraint: rl[i] <= aij * x_j <= ru[i]
        Real impl_lb, impl_ub;
        if (aij > 0.0) {
            impl_lb = is_inf(-rl[i]) ? -INF : rl[i] / aij;
            impl_ub = is_inf( ru[i]) ?  INF : ru[i] / aij;
        } else {
            impl_lb = is_inf( ru[i]) ? -INF : ru[i] / aij;
            impl_ub = is_inf(-rl[i]) ?  INF : rl[i] / aij;
        }

        Real tight_lb = std::max(cl[j], impl_lb);
        Real tight_ub = std::min(cu[j], impl_ub);

        if (tight_lb > tight_ub + EPS_FEASIBILITY) {
            infeasible = true;
            return true;
        }

        PresolveRecord rec;
        rec.rule     = PresolveRuleType::SINGLETON_ROW;
        rec.row_idx  = i;
        rec.col_idx  = j;
        rec.a_ij     = aij;
        rec.saved_lb = cl[j];
        rec.saved_ub = cu[j];
        rec.saved_rl = rl[i];
        rec.saved_ru = ru[i];
        s.push(rec);

        new_cl[j] = tight_lb;
        new_cu[j] = tight_ub;
        row_active[i] = false;
        changed = true;
    }

    if (changed) {
        Problem tmp = p;
        tmp.set_col_bounds(std::move(new_cl), std::move(new_cu));
        std::vector<bool>  col_active(p.num_cols(), true);
        std::vector<Index> m1, m2;
        p = compact(tmp, row_active, col_active, m1, m2);
    }
    return changed;
}

// ── Reduction: Singleton columns ──────────────────────────────────────────────

bool Presolve::reduce_singleton_cols(Problem& p, PresolveStack& s) {
    const auto& A  = p.constraint_matrix();
    const auto& cl = p.col_lower();
    const auto& cu = p.col_upper();
    const auto& c  = p.objective();
    const auto& rl = p.row_lower();
    const auto& ru = p.row_upper();

    bool changed = false;
    std::vector<bool> col_active(p.num_cols(), true);
    std::vector<bool> row_active(p.num_rows(), true);
    std::vector<Real> new_rl = rl, new_ru = ru;
    Real obj_offset = 0.0;

    for (Index j = 0; j < p.num_cols(); ++j) {
        if (A.col_length(j) != 1) continue;

        Index k   = A.col_begin(j);
        Index i   = A.row_index[k];
        Real  aij = A.values[k];

        if (!row_active[i]) continue;
        if (std::abs(aij) < EPS_ZERO) continue;

        // Determine optimal x_j value
        Real cj     = c[j];
        Real x_opt;

        if (cj > EPS_ZERO) {
            // Minimizing and c_j > 0: want x_j as small as possible
            x_opt = cl[j];
        } else if (cj < -EPS_ZERO) {
            // c_j < 0: want x_j as large as possible
            if (is_inf(cu[j])) return changed; // skip — would be unbounded
            x_opt = cu[j];
        } else {
            // c_j ≈ 0: choose to maximise slack in the constraint
            // x_j feasible range from constraint:
            Real con_lb = is_inf(-new_rl[i]) ? -INF : new_rl[i] / aij;
            Real con_ub = is_inf( new_ru[i]) ?  INF : new_ru[i] / aij;
            if (aij < 0) std::swap(con_lb, con_ub);
            x_opt = std::max(cl[j], con_lb);
            x_opt = std::min(x_opt, cu[j]);
            x_opt = std::min(x_opt, con_ub);
        }

        PresolveRecord rec;
        rec.rule      = PresolveRuleType::SINGLETON_COL;
        rec.row_idx   = i;
        rec.col_idx   = j;
        rec.a_ij      = aij;
        rec.c_j       = cj;
        rec.fix_value = x_opt;
        rec.saved_rl  = new_rl[i];
        rec.saved_ru  = new_ru[i];
        rec.saved_lb  = cl[j];
        rec.saved_ub  = cu[j];
        s.push(rec);

        // Substitute into row bounds
        Real sub = aij * x_opt;
        if (!is_inf(-new_rl[i])) new_rl[i] -= sub;
        if (!is_inf( new_ru[i])) new_ru[i] -= sub;

        obj_offset   += cj * x_opt;
        col_active[j] = false;
        changed = true;
    }

    if (changed) {
        Problem tmp = p;
        tmp.set_row_bounds(std::move(new_rl), std::move(new_ru));
        std::vector<Index> m1, m2;
        p = compact(tmp, row_active, col_active, m1, m2);
    }
    return changed;
}

// ── Reduction: Forcing rows ───────────────────────────────────────────────────

bool Presolve::reduce_forcing_rows(Problem& p, PresolveStack& s, bool& infeasible) {
    std::vector<Real> min_act, max_act;
    compute_row_activities(p, min_act, max_act);

    const auto& A  = p.constraint_matrix();
    const auto& cl = p.col_lower();
    const auto& cu = p.col_upper();
    const auto& rl = p.row_lower();
    const auto& ru = p.row_upper();

    bool changed = false;
    std::vector<bool> row_active(p.num_rows(), true);
    std::vector<Real> new_cl = cl, new_cu = cu;

    for (Index i = 0; i < p.num_rows(); ++i) {
        // Infeasibility check
        if (min_act[i] > ru[i] + EPS_FEASIBILITY ||
            max_act[i] < rl[i] - EPS_FEASIBILITY) {
            infeasible = true;
            return true;
        }

        // Forcing to upper bound (max_activity <= row_upper)
        if (max_act[i] <= ru[i] + EPS_FEASIBILITY) {
            PresolveRecord rec;
            rec.rule    = PresolveRuleType::FORCING_ROW;
            rec.row_idx = i;

            for (Index k = A.col_begin(i); k < A.col_end(i); ++k) {
                Index j   = A.row_index[k];
                Real  aij = A.values[k];
                Real  fix = (aij > 0) ? cu[j] : cl[j];
                rec.affected_cols.push_back(j);
                rec.saved_col_lbs.push_back(cl[j]);
                rec.saved_col_ubs.push_back(cu[j]);
                new_cl[j] = fix;
                new_cu[j] = fix;
            }
            s.push(rec);
            row_active[i] = false;
            changed = true;
            continue;
        }

        // Forcing to lower bound (min_activity >= row_lower)
        if (min_act[i] >= rl[i] - EPS_FEASIBILITY) {
            PresolveRecord rec;
            rec.rule    = PresolveRuleType::FORCING_ROW;
            rec.row_idx = i;

            for (Index k = A.col_begin(i); k < A.col_end(i); ++k) {
                Index j   = A.row_index[k];
                Real  aij = A.values[k];
                Real  fix = (aij > 0) ? cl[j] : cu[j];
                rec.affected_cols.push_back(j);
                rec.saved_col_lbs.push_back(cl[j]);
                rec.saved_col_ubs.push_back(cu[j]);
                new_cl[j] = fix;
                new_cu[j] = fix;
            }
            s.push(rec);
            row_active[i] = false;
            changed = true;
        }
    }

    if (changed) {
        Problem tmp = p;
        tmp.set_col_bounds(std::move(new_cl), std::move(new_cu));
        std::vector<bool>  col_active(p.num_cols(), true);
        std::vector<Index> m1, m2;
        p = compact(tmp, row_active, col_active, m1, m2);
    }
    return changed;
}

// ── Reduction: Implied bound tightening ──────────────────────────────────────

bool Presolve::tighten_implied_bounds(Problem& p, PresolveStack& s, bool& infeasible) {
    std::vector<Real> min_act, max_act;
    compute_row_activities(p, min_act, max_act);

    const auto& A  = p.constraint_matrix();
    const auto& cl = p.col_lower();
    const auto& cu = p.col_upper();
    const auto& rl = p.row_lower();
    const auto& ru = p.row_upper();

    bool changed = false;
    std::vector<Real> new_cl = cl, new_cu = cu;

    for (Index j = 0; j < p.num_cols(); ++j) {
        for (Index k = A.col_begin(j); k < A.col_end(j); ++k) {
            Index i   = A.row_index[k];
            Real  aij = A.values[k];
            if (std::abs(aij) < EPS_ZERO) continue;

            // Contribution of x_j to row i
            Real contrib_min = (aij > 0) ? aij * cl[j] : aij * cu[j];
            Real contrib_max = (aij > 0) ? aij * cu[j] : aij * cl[j];

            // Implied bounds from upper constraint:  aij * x_j <= ru[i] - (max_act[i] - contrib_max)
            if (!is_inf(ru[i])) {
                Real slack = ru[i] - (max_act[i] - contrib_max);
                // aij * x_j <= slack
                if (aij > 0) {
                    Real impl_ub = slack / aij;
                    if (impl_ub < new_cu[j] - EPS_ZERO) {
                        if (impl_ub < new_cl[j] - EPS_FEASIBILITY) {
                            infeasible = true; return true;
                        }
                        PresolveRecord rec;
                        rec.rule    = PresolveRuleType::IMPLIED_BOUND;
                        rec.col_idx = j;
                        rec.saved_ub = new_cu[j];
                        s.push(rec);
                        new_cu[j] = impl_ub;
                        changed = true;
                    }
                } else {
                    Real impl_lb = slack / aij;   // aij < 0 flips direction
                    if (impl_lb > new_cl[j] + EPS_ZERO) {
                        if (impl_lb > new_cu[j] + EPS_FEASIBILITY) {
                            infeasible = true; return true;
                        }
                        PresolveRecord rec;
                        rec.rule    = PresolveRuleType::IMPLIED_BOUND;
                        rec.col_idx = j;
                        rec.saved_lb = new_cl[j];
                        s.push(rec);
                        new_cl[j] = impl_lb;
                        changed = true;
                    }
                }
            }

            // Implied bounds from lower constraint: aij * x_j >= rl[i] - (min_act[i] - contrib_min)
            if (!is_inf(-rl[i])) {
                Real slack = rl[i] - (min_act[i] - contrib_min);
                if (aij > 0) {
                    Real impl_lb = slack / aij;
                    if (impl_lb > new_cl[j] + EPS_ZERO) {
                        if (impl_lb > new_cu[j] + EPS_FEASIBILITY) {
                            infeasible = true; return true;
                        }
                        PresolveRecord rec;
                        rec.rule    = PresolveRuleType::IMPLIED_BOUND;
                        rec.col_idx = j;
                        rec.saved_lb = new_cl[j];
                        s.push(rec);
                        new_cl[j] = impl_lb;
                        changed = true;
                    }
                } else {
                    Real impl_ub = slack / aij;
                    if (impl_ub < new_cu[j] - EPS_ZERO) {
                        if (impl_ub < new_cl[j] - EPS_FEASIBILITY) {
                            infeasible = true; return true;
                        }
                        PresolveRecord rec;
                        rec.rule    = PresolveRuleType::IMPLIED_BOUND;
                        rec.col_idx = j;
                        rec.saved_ub = new_cu[j];
                        s.push(rec);
                        new_cu[j] = impl_ub;
                        changed = true;
                    }
                }
            }
        }
    }

    if (changed)
        p.set_col_bounds(std::move(new_cl), std::move(new_cu));

    return changed;
}

// ── Reduction: Duplicate rows ─────────────────────────────────────────────────

bool Presolve::reduce_duplicate_rows(Problem& p, PresolveStack& s) {
    // Two rows are duplicates if their sparsity patterns and coefficient
    // ratios are identical. Merge by tightening bounds.
    // Simple O(m^2) check — adequate for typical presolve sizes.
    const auto& A  = p.constraint_matrix();
    const auto& rl = p.row_lower();
    const auto& ru = p.row_upper();

    bool changed = false;
    std::vector<bool> row_active(p.num_rows(), true);
    std::vector<Real> new_rl = rl, new_ru = ru;

    for (Index i = 0; i < p.num_rows(); ++i) {
        if (!row_active[i]) continue;
        for (Index ii = i + 1; ii < p.num_rows(); ++ii) {
            if (!row_active[ii]) continue;
            // Quick length check
            // (full duplicate detection deferred to advanced version)
        }
    }
    return changed; // stub — will be fully expanded
}

bool Presolve::reduce_dominated_cols(Problem& p, PresolveStack& s) {
    return false; // stub
}

bool Presolve::reduce_duplicate_cols(Problem& p, PresolveStack& s) {
    return false; // stub
}

// ── Reduction: Coefficient tightening (MILP) ─────────────────────────────────

bool Presolve::tighten_coefficients(Problem& p, PresolveStack& s) {
    const auto& A   = p.constraint_matrix();
    const auto& rl  = p.row_lower();
    const auto& ru  = p.row_upper();
    const auto& vt  = p.var_types();

    bool changed = false;
    // Access mutable versions through copies
    std::vector<Real> new_rl = rl, new_ru = ru;

    for (Index i = 0; i < p.num_rows(); ++i) {
        for (Index k = A.col_begin(i); k < A.col_end(i); ++k) {
            Index j   = A.row_index[k];
            Real  aij = A.values[k];
            if (vt[j] != VarType::INTEGER && vt[j] != VarType::BINARY) continue;
            if (std::abs(aij) < EPS_ZERO) continue;

            // Example: a_ij * x_j <= rhs, x_j integer
            // floor(rhs / a_ij) gives tighter integer bound
            if (!is_inf(ru[i]) && aij > 0) {
                Real rhs_per_unit = ru[i] / aij;
                Real floor_bound  = std::floor(rhs_per_unit + EPS_FEASIBILITY);
                if (floor_bound * aij < ru[i] - EPS_FEASIBILITY) {
                    PresolveRecord rec;
                    rec.rule     = PresolveRuleType::COEFFICIENT_TIGHTENING;
                    rec.row_idx  = i;
                    rec.saved_ru = new_ru[i];
                    s.push(rec);
                    new_ru[i] = floor_bound * aij;
                    changed = true;
                }
            }
        }
    }

    if (changed)
        p.set_row_bounds(std::move(new_rl), std::move(new_ru));

    return changed;
}

// ── Reduction: Probing (binary variables) ────────────────────────────────────

bool Presolve::run_probing(Problem& p, PresolveStack& s, bool& infeasible) {
    const auto& vt = p.var_types();
    const auto& cl = p.col_lower();
    const auto& cu = p.col_upper();

    bool changed = false;

    for (Index j = 0; j < p.num_cols(); ++j) {
        if (vt[j] != VarType::BINARY &&
            !(vt[j] == VarType::INTEGER &&
              nearly_equal(cl[j], 0.0) && nearly_equal(cu[j], 1.0))) continue;

        // Save bounds
        std::vector<Real> saved_cl = p.col_lower();
        std::vector<Real> saved_cu = p.col_upper();

        // ── Probe x_j = 0 ──
        {
            std::vector<Real> probe_cl = saved_cl;
            std::vector<Real> probe_cu = saved_cu;
            probe_cu[j] = 0.0;
            Problem p0 = p;
            p0.set_col_bounds(std::move(probe_cl), std::move(probe_cu));
            bool inf0 = false;
            propagate_bounds(p0, inf0);

            // ── Probe x_j = 1 ──
            std::vector<Real> probe_cl1 = saved_cl;
            std::vector<Real> probe_cu1 = saved_cu;
            probe_cl1[j] = 1.0;
            Problem p1 = p;
            p1.set_col_bounds(std::move(probe_cl1), std::move(probe_cu1));
            bool inf1 = false;
            propagate_bounds(p1, inf1);

            if (inf0 && inf1) { infeasible = true; return true; }

            if (inf0) {
                // x_j must be 1
                PresolveRecord rec;
                rec.rule      = PresolveRuleType::PROBING;
                rec.col_idx   = j;
                rec.fix_value = 1.0;
                rec.saved_lb  = cl[j];
                rec.saved_ub  = cu[j];
                s.push(rec);
                std::vector<Real> new_cl = p.col_lower();
                std::vector<Real> new_cu = p.col_upper();
                new_cl[j] = 1.0;
                p.set_col_bounds(std::move(new_cl), std::move(new_cu));
                changed = true;
            } else if (inf1) {
                // x_j must be 0
                PresolveRecord rec;
                rec.rule      = PresolveRuleType::PROBING;
                rec.col_idx   = j;
                rec.fix_value = 0.0;
                rec.saved_lb  = cl[j];
                rec.saved_ub  = cu[j];
                s.push(rec);
                std::vector<Real> new_cl = p.col_lower();
                std::vector<Real> new_cu = p.col_upper();
                new_cu[j] = 0.0;
                p.set_col_bounds(std::move(new_cl), std::move(new_cu));
                changed = true;
            } else {
                // Tighten bounds agreed by both probings
                const auto& cl0 = p0.col_lower(); const auto& cu0 = p0.col_upper();
                const auto& cl1 = p1.col_lower(); const auto& cu1 = p1.col_upper();
                std::vector<Real> new_cl = p.col_lower();
                std::vector<Real> new_cu = p.col_upper();
                bool local_changed = false;
                for (Index k = 0; k < p.num_cols(); ++k) {
                    Real agreed_lb = std::max(cl0[k], cl1[k]);
                    Real agreed_ub = std::min(cu0[k], cu1[k]);
                    if (agreed_lb > new_cl[k] + EPS_FEASIBILITY) {
                        new_cl[k] = agreed_lb; local_changed = true;
                    }
                    if (agreed_ub < new_cu[k] - EPS_FEASIBILITY) {
                        new_cu[k] = agreed_ub; local_changed = true;
                    }
                }
                if (local_changed) {
                    p.set_col_bounds(std::move(new_cl), std::move(new_cu));
                    changed = true;
                }
            }
        }
    }
    return changed;
}

// ── Internal helpers ──────────────────────────────────────────────────────────

void Presolve::compute_row_activities(const Problem& p,
                                       std::vector<Real>& min_act,
                                       std::vector<Real>& max_act) const {
    const auto& A  = p.constraint_matrix();
    const auto& cl = p.col_lower();
    const auto& cu = p.col_upper();

    min_act.assign(p.num_rows(), 0.0);
    max_act.assign(p.num_rows(), 0.0);

    for (Index j = 0; j < p.num_cols(); ++j) {
        for (Index k = A.col_begin(j); k < A.col_end(j); ++k) {
            Index i   = A.row_index[k];
            Real  aij = A.values[k];
            if (aij > 0) {
                if (!is_inf(-cl[j])) min_act[i] += aij * cl[j]; else min_act[i] = -INF;
                if (!is_inf( cu[j])) max_act[i] += aij * cu[j]; else max_act[i] =  INF;
            } else {
                if (!is_inf( cu[j])) min_act[i] += aij * cu[j]; else min_act[i] = -INF;
                if (!is_inf(-cl[j])) max_act[i] += aij * cl[j]; else max_act[i] =  INF;
            }
        }
    }
}

bool Presolve::propagate_bounds(Problem& p, bool& infeasible) {
    // Single pass of implied bound tightening (used by probing)
    PresolveStack dummy;
    return tighten_implied_bounds(p, dummy, infeasible);
}

Problem Presolve::compact(const Problem& p,
                           const std::vector<bool>& row_active,
                           const std::vector<bool>& col_active,
                           std::vector<Index>& old_to_new_row,
                           std::vector<Index>& old_to_new_col) const {
    // Build index maps
    old_to_new_row.assign(p.num_rows(), -1);
    old_to_new_col.assign(p.num_cols(), -1);
    Index new_m = 0, new_n = 0;
    for (Index i = 0; i < p.num_rows(); ++i)
        if (row_active[i]) old_to_new_row[i] = new_m++;
    for (Index j = 0; j < p.num_cols(); ++j)
        if (col_active[j]) old_to_new_col[j] = new_n++;

    // Build new matrix in triplet form
    const auto& A = p.constraint_matrix();
    std::vector<Index> tri_r, tri_c;
    std::vector<Real>  tri_v;

    for (Index j = 0; j < p.num_cols(); ++j) {
        if (!col_active[j]) continue;
        Index nj = old_to_new_col[j];
        for (Index k = A.col_begin(j); k < A.col_end(j); ++k) {
            Index i = A.row_index[k];
            if (!row_active[i]) continue;
            tri_r.push_back(old_to_new_row[i]);
            tri_c.push_back(nj);
            tri_v.push_back(A.values[k]);
        }
    }

    Problem q;
    q.set_dimensions(new_m, new_n);

    // Objective, types, col bounds
    const auto& c  = p.objective();
    const auto& cl = p.col_lower();
    const auto& cu = p.col_upper();
    const auto& vt = p.var_types();
    std::vector<Real>    new_c(new_n), new_cl(new_n), new_cu(new_n);
    std::vector<VarType> new_vt(new_n);
    for (Index j = 0; j < p.num_cols(); ++j) {
        if (!col_active[j]) continue;
        Index nj   = old_to_new_col[j];
        new_c[nj]  = c[j];
        new_cl[nj] = cl[j];
        new_cu[nj] = cu[j];
        new_vt[nj] = vt[j];
    }
    q.set_objective(std::move(new_c));
    q.set_col_bounds(std::move(new_cl), std::move(new_cu));
    q.set_var_types(std::move(new_vt));

    // Row bounds
    const auto& rl = p.row_lower();
    const auto& ru = p.row_upper();
    std::vector<Real> new_rl(new_m), new_ru(new_m);
    for (Index i = 0; i < p.num_rows(); ++i) {
        if (!row_active[i]) continue;
        Index ni    = old_to_new_row[i];
        new_rl[ni]  = rl[i];
        new_ru[ni]  = ru[i];
    }
    q.set_row_bounds(std::move(new_rl), std::move(new_ru));
    q.set_objective_sense(p.sense());

    // Constraint matrix
    auto new_A = SparseMatrixCSC::from_triplets(new_m, new_n, tri_r, tri_c, tri_v);
    q.set_constraint_matrix(std::move(new_A));

    return q;
}

} // namespace suplex
