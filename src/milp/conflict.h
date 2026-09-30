#pragma once
#include "../core/types.h"
#include "../core/problem.h"
#include "../core/solution.h"
#include "node.h"
#include "cuts.h"
#include <vector>

namespace suplex {

// ── Conflict analysis ─────────────────────────────────────────────────────────
///
/// When a B&B node's LP is infeasible, analyse WHY using the Farkas ray
/// and generate a "conflict cut" (no-good clause) to prune similar subtrees.
///
class ConflictAnalysis {
public:
    /// Analyse the infeasible node and produce a conflict cut.
    ///
    /// @param farkas_ray  Farkas dual certificate (y^T A >= 0, y^T b < 0)
    /// @param node        The infeasible node
    /// @param branch_path Sequence of (col, bound_change) from root to node
    /// @returns           A conflict cut, or empty if analysis fails
    std::optional<Cut> analyse(
        const std::vector<Real>&  farkas_ray,
        const BBNode&             node,
        const std::vector<std::pair<Index, BoundChange>>& branch_path,
        const Problem&            p);

private:
    /// Identify which columns contribute to the Farkas infeasibility.
    std::vector<Index> responsible_cols(
        const std::vector<Real>&  farkas_ray,
        const Problem&            p) const;

    /// Remove redundant decisions from the conflict set.
    std::vector<std::pair<Index, BoundChange>> minimize_conflict(
        const std::vector<std::pair<Index, BoundChange>>& conflict_set,
        const std::vector<Index>& responsible_cols) const;

    /// Convert a minimal conflict set into a no-good cut.
    Cut build_nogood_cut(
        const std::vector<std::pair<Index, BoundChange>>& minimal_conflict,
        const Problem& p) const;
};

} // namespace suplex
