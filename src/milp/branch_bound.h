#pragma once
#include "../core/types.h"
#include "../core/problem.h"
#include "../core/solution.h"
#include "../simplex/lp_solver.h"
#include "node.h"
#include "presolve.h"
#include <vector>
#include <memory>
#include <functional>

namespace suplex {

// ── MILP callback ─────────────────────────────────────────────────────────────
/// Override to inspect or control the B&B solve.
class MILPCallback {
public:
    virtual ~MILPCallback() = default;
    virtual void on_new_incumbent(const Solution& sol) {}
    virtual void on_node_solved(int64_t node_id, Real bound) {}
    virtual bool should_terminate() { return false; }
};

// ── Node selection strategy ───────────────────────────────────────────────────
enum class NodeSelection : uint8_t {
    BEST_FIRST,
    DEPTH_FIRST,
    BEST_ESTIMATE,
    HYBRID        ///< Default: plunge depth-first, then best-first
};

// ── Branch-and-Bound tree manager ────────────────────────────────────────────
///
/// Orchestrates the complete MILP solve:
///   root LP → cuts → heuristics → B&B loop → postsolve
///
class BranchBound {
public:
    explicit BranchBound(std::shared_ptr<LPSolver> lp_solver);

    // ── Main solve entry point ────────────────────────────────────────────
    SolverStatus solve(const Problem& problem, Solution& solution);

    // ── Configuration ─────────────────────────────────────────────────────
    void set_gap_tolerance(Real gap)              { gap_tol_ = gap; }
    void set_node_limit(int64_t limit)            { node_limit_ = limit; }
    void set_time_limit(Real seconds)             { time_limit_ = seconds; }
    void set_log_level(LogLevel level)            { log_level_ = level; }
    void set_node_selection(NodeSelection s)      { node_selection_ = s; }
    void set_callback(std::shared_ptr<MILPCallback> cb) { callback_ = cb; }
    void set_cut_rounds_at_root(int n)            { cut_rounds_root_ = n; }
    void set_max_cut_rounds_per_node(int n)       { cut_rounds_node_ = n; }

    const TreeStats& tree_stats() const { return stats_; }

private:
    std::shared_ptr<LPSolver>     lp_solver_;
    std::shared_ptr<MILPCallback> callback_;

    // Configuration
    Real          gap_tol_           = 1e-4;
    int64_t       node_limit_        = 1000000000LL;
    Real          time_limit_        = 3600.0;
    LogLevel      log_level_         = LogLevel::INFO;
    NodeSelection node_selection_    = NodeSelection::HYBRID;
    int           cut_rounds_root_   = 20;
    int           cut_rounds_node_   = 5;

    // Runtime state
    TreeStats     stats_;
    Real          incumbent_obj_     = 1e30;    ///< Best integer solution found
    Solution      incumbent_sol_;
    double        start_time_        = 0.0;

    // ── Internal phases ───────────────────────────────────────────────────
    /// Solve root node LP and return status.
    SolverStatus solve_root_lp(const Problem& p, Solution& sol);

    /// Separate cuts at a node. Returns number of cuts added.
    int separate_cuts(const Problem& p, const Solution& lp_sol,
                      Problem& p_with_cuts);

    /// Run primal heuristics. Returns true if new incumbent found.
    bool run_heuristics(const Problem& p, const Solution& lp_sol);

    /// Process one B&B node: solve LP, prune/branch.
    void process_node(BBNode& node, const Problem& base_problem,
                      std::vector<std::unique_ptr<BBNode>>& new_nodes);

    /// Select next node from pool according to node_selection_ strategy.
    BBNode* select_node(std::vector<std::unique_ptr<BBNode>>& pool,
                        const BBNode* last_node);

    /// Check if all integer variables are integral in `sol`.
    bool is_integer_feasible(const Problem& p, const Solution& sol) const;

    /// Update incumbent if sol is better.
    void try_update_incumbent(const Problem& p, const Solution& sol,
                              Index base_rows);

    /// Apply bound changes from a B&B node to a problem copy.
    void apply_bound_changes(Problem& p,
                             const std::vector<BoundChange>& changes) const;

    /// Compute current MIP gap.
    Real compute_gap(Real best_bound) const;

    /// Log one-line progress update.
    void log_progress(int64_t node_id, Real lp_bound,
                      int depth, int open_nodes) const;

    void remove_from_pool(std::vector<std::unique_ptr<BBNode>>& pool,
                          BBNode* node);
};

} // namespace suplex
