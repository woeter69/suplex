#include "conflict.h"
#include <algorithm>
#include <cmath>

namespace suplex {

std::optional<Cut> ConflictAnalysis::analyse(
        const std::vector<Real>&  farkas_ray,
        const BBNode&             node,
        const std::vector<std::pair<Index, BoundChange>>& branch_path,
        const Problem&            p) {

    if (farkas_ray.empty()) return std::nullopt;

    // Step 1: identify which columns appear in responsible rows
    auto resp_cols = responsible_cols(farkas_ray, p);
    if (resp_cols.empty()) return std::nullopt;

    // Step 2: filter branch_path to decisions that affect responsible columns
    std::vector<std::pair<Index, BoundChange>> conflict_set;
    for (const auto& [col, bc] : branch_path) {
        bool relevant = std::find(resp_cols.begin(), resp_cols.end(), col)
                        != resp_cols.end();
        if (relevant) conflict_set.push_back({col, bc});
    }
    if (conflict_set.empty()) return std::nullopt;

    // Step 3: minimise (remove redundant decisions)
    auto minimal = minimize_conflict(conflict_set, resp_cols);
    if (minimal.empty()) return std::nullopt;

    // Step 4: build no-good cut
    return build_nogood_cut(minimal, p);
}

std::vector<Index> ConflictAnalysis::responsible_cols(
        const std::vector<Real>& farkas_ray,
        const Problem&           p) const {
    // Identify rows with nonzero Farkas multiplier
    std::vector<Index> resp_rows;
    for (Index i = 0; i < static_cast<Index>(farkas_ray.size()); ++i)
        if (std::abs(farkas_ray[i]) > EPS_ZERO) resp_rows.push_back(i);

    // Collect all columns that appear in those rows
    const auto& A = p.constraint_matrix();
    std::vector<Index> resp_cols;
    for (Index j = 0; j < p.num_cols(); ++j)
        for (Index k = A.col_begin(j); k < A.col_end(j); ++k)
            if (std::find(resp_rows.begin(), resp_rows.end(), A.row_index[k])
                    != resp_rows.end()) {
                resp_cols.push_back(j);
                break;
            }

    return resp_cols;
}

std::vector<std::pair<Index, BoundChange>> ConflictAnalysis::minimize_conflict(
        const std::vector<std::pair<Index, BoundChange>>& conflict_set,
        const std::vector<Index>& /*responsible_cols*/) const {
    // Conservative: return all decisions (full minimisation is expensive)
    // In production this would iteratively remove redundant decisions
    return conflict_set;
}

Cut ConflictAnalysis::build_nogood_cut(
        const std::vector<std::pair<Index, BoundChange>>& minimal_conflict,
        const Problem& p) const {
    // Build a no-good clause:
    //   Σ (x_j if branched down) + Σ (1 - x_j if branched up) >= 1
    // For binary variables: prevents revisiting the same combination of fixings
    Cut cut;
    cut.type = CutType::GOMORY;   // re-use as a general cut type
    cut.rhs  = 1.0;
    cut.name = "conflict";

    const auto& cu = p.col_upper();
    const auto& cl = p.col_lower();

    for (const auto& [col, bc] : minimal_conflict) {
        // If branched DOWN (x_j <= floor): literal is x_j
        // If branched UP   (x_j >= ceil) : literal is 1 - x_j  → coefficient -1, adjust rhs
        bool down_branch = (bc.new_upper < cu[col] - EPS_ZERO);
        if (down_branch) {
            // x_j <= floor: add (x_j - lb) / (ub - lb) to handle general integers
            cut.indices.push_back(col);
            cut.coefficients.push_back(1.0);
        } else {
            // x_j >= ceil: add (ub - x_j) / (ub - lb) → coefficient -1, adjust rhs
            cut.indices.push_back(col);
            cut.coefficients.push_back(-1.0);
            cut.rhs -= 1.0;   // accounts for -x_j term
        }
    }
    return cut;
}

} // namespace suplex
