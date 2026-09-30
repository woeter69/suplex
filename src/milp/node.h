#pragma once
#include "../core/types.h"
#include <vector>
#include <cstdint>

namespace suplex {

// ── B&B node ──────────────────────────────────────────────────────────────────

enum class BranchDirection : uint8_t { DOWN, UP };

enum class NodeStatus : uint8_t {
    OPEN,
    SOLVING,
    SOLVED,
    PRUNED_BOUND,
    PRUNED_INFEASIBLE,
    BRANCHED,
    INTEGER_FEASIBLE
};

/// One bound change applied at this node (difference from parent).
struct BoundChange {
    Index col_idx;
    Real  new_lower;
    Real  new_upper;
};

/// Single node in the branch-and-bound search tree.
/// Stores only DELTAS from the parent — not a full problem copy.
struct BBNode {
    int64_t  id        = 0;
    int64_t  parent_id = -1;     ///< -1 for root
    int      depth     = 0;

    // Branching decision that created this node
    Index          branch_var   = -1;
    Real           branch_value = 0.0;  ///< Fractional value that was branched on
    BranchDirection branch_dir  = BranchDirection::DOWN;

    // Bound changes from parent (usually just 1 — the branching bound)
    std::vector<BoundChange> bound_changes;

    // LP relaxation result at this node
    Real         lp_bound  = -1e30;
    NodeStatus   status    = NodeStatus::OPEN;

    // Warm-start basis (row-indexed: basic_vars[i] = column of i-th basic var)
    std::vector<Index> basis;

    // Estimated integer objective (used by best-estimate node selection)
    Real estimate = -1e30;

    // ── Comparators for priority queues ──────────────────────────────────────
    /// Best-first: lowest bound for minimization.
    struct CmpBestBound {
        bool operator()(const BBNode* a, const BBNode* b) const {
            return a->lp_bound > b->lp_bound;  // min-heap
        }
    };
    /// Depth-first: deepest node first.
    struct CmpDepthFirst {
        bool operator()(const BBNode* a, const BBNode* b) const {
            return a->depth < b->depth;
        }
    };
    /// Best-estimate: uses pseudo-cost estimate of integer objective.
    struct CmpBestEstimate {
        bool operator()(const BBNode* a, const BBNode* b) const {
            return a->estimate > b->estimate;  // min-heap on estimate
        }
    };
};

// ── Tree statistics (logged periodically) ────────────────────────────────────
struct TreeStats {
    int64_t nodes_created    = 0;
    int64_t nodes_solved     = 0;
    int64_t nodes_pruned     = 0;
    int64_t nodes_infeasible = 0;
    int64_t cuts_added       = 0;
    int64_t heuristic_sols   = 0;
    int     max_depth        = 0;
    double  nodes_per_second = 0.0;
};

} // namespace suplex
