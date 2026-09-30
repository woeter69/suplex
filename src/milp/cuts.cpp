#include "cuts.h"
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>

namespace suplex {

// ── CutPool ───────────────────────────────────────────────────────────────────

void CutPool::add(Cut cut) { cuts_.push_back(std::move(cut)); }

void CutPool::age_all() { for (auto& c : cuts_) c.age++; }

void CutPool::remove_old(int max_age) {
    cuts_.erase(std::remove_if(cuts_.begin(), cuts_.end(),
                               [max_age](const Cut& c) { return c.age > max_age; }),
                cuts_.end());
}

// ── CutGenerator ─────────────────────────────────────────────────────────────

CutGenerator::CutGenerator(Index num_cols) : num_cols_(num_cols) {}

// ── Efficacy: normalised distance from LP solution to cut hyperplane ──────────

Real CutGenerator::compute_efficacy(const Cut& cut, const Solution& lp_sol) const {
    Real lhs = 0.0;
    for (Index k = 0; k < static_cast<Index>(cut.indices.size()); ++k)
        lhs += cut.coefficients[k] * lp_sol.primal_values[cut.indices[k]];

    Real violation = cut.rhs - lhs;  // >= 0 if violated (we want >= rhs)
    Real norm = 0.0;
    for (Real v : cut.coefficients) norm += v * v;
    norm = std::sqrt(norm);
    return (norm > EPS_ZERO) ? violation / norm : 0.0;
}

Real CutGenerator::compute_orthogonality(const Cut& a, const Cut& b) const {
    // Cosine between cut normals — high = parallel (bad), low = orthogonal (good)
    Real dot = 0.0, na = 0.0, nb = 0.0;

    // Build sparse dot product
    // O(|a| * |b|) — fine for typical small cuts
    std::vector<Real> dense_b(num_cols_, 0.0);
    for (Index k = 0; k < static_cast<Index>(b.indices.size()); ++k)
        dense_b[b.indices[k]] = b.coefficients[k];

    for (Index k = 0; k < static_cast<Index>(a.indices.size()); ++k) {
        dot += a.coefficients[k] * dense_b[a.indices[k]];
        na  += a.coefficients[k] * a.coefficients[k];
    }
    for (Real v : b.coefficients) nb += v * v;

    if (na < EPS_ZERO || nb < EPS_ZERO) return 1.0;
    Real cosine = dot / (std::sqrt(na) * std::sqrt(nb));
    return std::abs(cosine);   // 0 = orthogonal (good), 1 = parallel (bad)
}

// ── Gomory cuts ───────────────────────────────────────────────────────────────

bool CutGenerator::extract_tableau_row(Index basic_pos,
                                        const std::vector<Index>& basis,
                                        const Problem& p,
                                        std::vector<Real>& row_coeffs,
                                        Real& rhs) const {
    // In a full implementation, this would call LPSolver::btran(e_{basic_pos})
    // and then dot with each column. For now we stub it.
    // Person 1 must provide: basis_inverse_row(basic_pos) -> dense vector
    // Then: tableau_coeff[j] = dot(basis_inv_row, A[:,j])
    (void)basic_pos; (void)basis; (void)p;
    row_coeffs.clear();
    rhs = 0.0;
    return false;   // stub — returns false until Person 1 integrates
}

std::vector<Cut> CutGenerator::gomory_cuts(const Problem&  p,
                                            const Solution& lp_sol,
                                            const std::vector<Index>& basis,
                                            int max_cuts) {
    std::vector<Cut> cuts;
    const auto& vt = p.var_types();

    for (Index pos = 0; pos < static_cast<Index>(basis.size()); ++pos) {
        Index j = basis[pos];
        if (j >= p.num_cols()) continue;               // skip slacks
        if (vt[j] == VarType::CONTINUOUS) continue;    // only integer vars

        Real val = lp_sol.primal_values[j];
        Real f0  = val - std::floor(val);

        if (f0 < 0.05 || f0 > 0.95) continue;  // near-integer — skip

        // Extract tableau row
        std::vector<Real> row_coeffs;
        Real row_rhs;
        if (!extract_tableau_row(pos, basis, p, row_coeffs, row_rhs)) continue;

        // Build Gomory cut: Σ gomory_coeff_j x_j >= 1
        Cut cut;
        cut.type = CutType::GOMORY;

        for (Index nj = 0; nj < p.num_cols(); ++nj) {
            if (std::abs(row_coeffs[nj]) < EPS_ZERO) continue;
            Real aij = row_coeffs[nj];
            Real fj  = aij - std::floor(aij);

            Real gc = (fj <= f0) ? fj / f0 : (1.0 - fj) / (1.0 - f0);
            if (std::abs(gc) < EPS_ZERO) continue;

            cut.indices.push_back(nj);
            cut.coefficients.push_back(gc);
        }
        cut.rhs = 1.0;

        Real eff = compute_efficacy(cut, lp_sol);
        if (eff < min_efficacy_) continue;
        cut.efficacy = eff;

        cuts.push_back(std::move(cut));
        if (static_cast<int>(cuts.size()) >= max_cuts) break;
    }
    return cuts;
}

// ── MIR cuts ──────────────────────────────────────────────────────────────────

Cut CutGenerator::mir_from_row(const std::vector<Real>& row_coeffs,
                                const std::vector<Index>& col_map,
                                Real rhs,
                                const Problem& p,
                                const Solution& lp_sol) const {
    Real f0 = rhs - std::floor(rhs);
    Cut cut;
    cut.type = CutType::MIR;

    if (f0 < 0.05 || f0 > 0.95) return cut;   // empty cut

    const auto& vt = p.var_types();

    for (Index k = 0; k < static_cast<Index>(col_map.size()); ++k) {
        Index j   = col_map[k];
        Real  aij = row_coeffs[k];
        if (std::abs(aij) < EPS_ZERO) continue;

        Real mir_coeff;
        if (vt[j] == VarType::CONTINUOUS) {
            // Continuous: coefficient ≥ 0 → aij/(1-f0), else 0
            mir_coeff = (aij > 0) ? aij / (1.0 - f0) : 0.0;
        } else {
            // Integer: MIR formula
            Real fj = aij - std::floor(aij);
            if (fj <= f0)
                mir_coeff = std::floor(aij);
            else
                mir_coeff = std::floor(aij) + (fj - f0) / (1.0 - f0);
        }

        if (std::abs(mir_coeff) > EPS_ZERO) {
            cut.indices.push_back(j);
            cut.coefficients.push_back(mir_coeff);
        }
    }
    cut.rhs = std::floor(rhs);
    return cut;
}

std::vector<Cut> CutGenerator::mir_cuts(const Problem&  p,
                                         const Solution& lp_sol,
                                         int max_cuts) {
    std::vector<Cut> cuts;
    const auto& A  = p.constraint_matrix();
    const auto& rl = p.row_lower();
    const auto& ru = p.row_upper();

    for (Index i = 0; i < p.num_rows(); ++i) {
        // Try the <= form of each constraint: Σ a_ij x_j <= ru[i]
        if (ru[i] >= INF) continue;

        std::vector<Real>  row_coeffs;
        std::vector<Index> col_map;

        for (Index k = A.col_begin(i); k < A.col_end(i); ++k) {
            col_map.push_back(A.row_index[k]);
            row_coeffs.push_back(A.values[k]);
        }

        Cut cut = mir_from_row(row_coeffs, col_map, ru[i], p, lp_sol);
        if (cut.indices.empty()) continue;

        Real eff = compute_efficacy(cut, lp_sol);
        if (eff < min_efficacy_) continue;
        cut.efficacy = eff;
        cuts.push_back(std::move(cut));
        if (static_cast<int>(cuts.size()) >= max_cuts) break;
    }
    return cuts;
}

// ── Knapsack cover cuts ───────────────────────────────────────────────────────

bool CutGenerator::find_minimal_cover(const std::vector<Index>& vars,
                                       const std::vector<Real>&  coeffs,
                                       Real                      rhs,
                                       std::vector<Index>&       cover) const {
    // Greedy: sort by coeff descending, add until sum > rhs
    std::vector<std::pair<Real, Index>> sorted;
    for (Index k = 0; k < static_cast<Index>(vars.size()); ++k)
        sorted.emplace_back(coeffs[k], k);
    std::sort(sorted.begin(), sorted.end(), std::greater<>());

    Real   sum = 0.0;
    cover.clear();
    for (auto& [coeff, k] : sorted) {
        cover.push_back(k);
        sum += coeff;
        if (sum > rhs + EPS_FEASIBILITY) return true;
    }
    return false;   // no cover found
}

Real CutGenerator::lifting_coefficient(const std::vector<Index>& cover,
                                        const std::vector<Real>&  coeffs,
                                        Index   new_var_idx, Real new_var_coeff,
                                        Real    knapsack_rhs) const {
    // Simplified sequential lifting: compute how much new_var can participate
    // without violating the cover inequality.
    Real cover_sum = 0.0;
    for (Index k : cover) cover_sum += coeffs[k];

    // Max alpha such that cover_cut + alpha*x_{new} <= |cover|-1 is valid
    // Conservative estimate: alpha = floor((cover_sum - knapsack_rhs) / new_var_coeff)
    if (new_var_coeff < EPS_ZERO) return 0.0;
    Real excess = cover_sum - knapsack_rhs;
    return std::max(0.0, std::floor(excess / new_var_coeff));
}

std::vector<Cut> CutGenerator::cover_cuts(const Problem&  p,
                                           const Solution& lp_sol,
                                           int max_cuts) {
    std::vector<Cut> cuts;
    const auto& A  = p.constraint_matrix();
    const auto& ru = p.row_upper();
    const auto& vt = p.var_types();
    const auto& cu = p.col_upper();

    for (Index i = 0; i < p.num_rows(); ++i) {
        if (ru[i] >= INF) continue;

        // Collect binary columns in this row
        std::vector<Index> bin_vars;
        std::vector<Real>  bin_coeffs;

        for (Index k = A.col_begin(i); k < A.col_end(i); ++k) {
            Index j = A.row_index[k];
            if (vt[j] != VarType::BINARY &&
                !(vt[j] == VarType::INTEGER && cu[j] <= 1.0 + EPS_ZERO)) continue;
            if (A.values[k] <= 0) continue;   // only positive coefficients
            bin_vars.push_back(j);
            bin_coeffs.push_back(A.values[k]);
        }

        if (bin_vars.size() < 2) continue;

        std::vector<Index> cover_idx;
        if (!find_minimal_cover(bin_vars, bin_coeffs, ru[i], cover_idx)) continue;

        // Build cover cut: Σ_{k in cover} x_{bin_vars[k]} <= |cover| - 1
        Cut cut;
        cut.type = CutType::COVER;
        cut.rhs  = static_cast<Real>(cover_idx.size()) - 1.0;

        for (Index k : cover_idx) {
            cut.indices.push_back(bin_vars[k]);
            cut.coefficients.push_back(1.0);
        }

        // Lift non-cover binary variables
        for (Index k = 0; k < static_cast<Index>(bin_vars.size()); ++k) {
            // Check if k is in cover
            bool in_cover = std::find(cover_idx.begin(), cover_idx.end(), k)
                            != cover_idx.end();
            if (in_cover) continue;

            Real alpha = lifting_coefficient(cover_idx, bin_coeffs,
                                             k, bin_coeffs[k], ru[i]);
            if (alpha > EPS_ZERO) {
                cut.indices.push_back(bin_vars[k]);
                cut.coefficients.push_back(alpha);
            }
        }

        Real eff = compute_efficacy(cut, lp_sol);
        if (eff < min_efficacy_) continue;
        cut.efficacy = eff;
        cuts.push_back(std::move(cut));
        if (static_cast<int>(cuts.size()) >= max_cuts) break;
    }
    return cuts;
}

// ── Clique cuts ───────────────────────────────────────────────────────────────

std::vector<Cut> CutGenerator::clique_cuts(const Problem&  p,
                                            const Solution& lp_sol,
                                            int max_cuts) {
    std::vector<Cut> cuts;
    const auto& A  = p.constraint_matrix();
    const auto& ru = p.row_upper();
    const auto& vt = p.var_types();
    const auto& cu = p.col_upper();

    // Simple clique detection: find rows where all binary vars sum to <= 1
    for (Index i = 0; i < p.num_rows(); ++i) {
        if (ru[i] > 1.0 + EPS_ZERO) continue;

        // All coefficients must be 1.0 and all vars binary
        std::vector<Index> clique;
        bool valid = true;
        for (Index k = A.col_begin(i); k < A.col_end(i); ++k) {
            Index j = A.row_index[k];
            if ((vt[j] != VarType::BINARY &&
                 !(vt[j] == VarType::INTEGER && cu[j] <= 1.0 + EPS_ZERO)) ||
                std::abs(A.values[k] - 1.0) > EPS_ZERO) {
                valid = false; break;
            }
            clique.push_back(j);
        }
        if (!valid || clique.size() < 2) continue;

        // The constraint itself IS the clique cut — check if violated
        Cut cut;
        cut.type = CutType::CLIQUE;
        cut.rhs  = 1.0;
        for (Index j : clique) {
            cut.indices.push_back(j);
            cut.coefficients.push_back(1.0);
        }

        Real eff = compute_efficacy(cut, lp_sol);
        if (eff < min_efficacy_) continue;
        cut.efficacy = eff;
        cuts.push_back(std::move(cut));
        if (static_cast<int>(cuts.size()) >= max_cuts) break;
    }
    return cuts;
}

// ── Cut selection ─────────────────────────────────────────────────────────────

std::vector<Cut> CutGenerator::select_cuts(std::vector<Cut>& candidates,
                                            int max_to_add) {
    // Score: efficacy - 0.3 * orthogonality_penalty
    // Sort by efficacy descending first
    std::sort(candidates.begin(), candidates.end(),
              [](const Cut& a, const Cut& b) { return a.efficacy > b.efficacy; });

    std::vector<Cut> selected;
    selected.reserve(max_to_add);

    for (auto& c : candidates) {
        if (static_cast<int>(selected.size()) >= max_to_add) break;
        if (c.efficacy < min_efficacy_) continue;

        // Check orthogonality against already selected cuts
        bool too_parallel = false;
        for (const auto& s : selected) {
            if (compute_orthogonality(c, s) > 0.9) {  // > 90% parallel → skip
                too_parallel = true;
                break;
            }
        }
        if (!too_parallel) selected.push_back(std::move(c));
    }
    return selected;
}

} // namespace suplex
