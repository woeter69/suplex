#pragma once
#include "../core/types.h"
#include "../core/problem.h"
#include "../core/solution.h"
#include <vector>
#include <memory>

namespace suplex {

// ── Presolve record types ────────────────────────────────────────────────────

/// Identifies which reduction rule was applied so postsolve can undo it.
enum class PresolveRuleType : uint8_t {
    EMPTY_ROW,
    SINGLETON_ROW,
    SINGLETON_COL,
    FORCING_ROW,
    DOMINATED_COL,
    DUPLICATE_ROW,
    DUPLICATE_COL,
    FREE_COL_SUBSTITUTION,
    IMPLIED_BOUND,
    FIXED_VARIABLE,
    PROBING,
    COEFFICIENT_TIGHTENING
};

/// One undo-able reduction. Stored on the PresolveStack.
struct PresolveRecord {
    PresolveRuleType rule;

    // Variable / constraint indices
    Index row_idx  = -1;
    Index col_idx  = -1;

    // Coefficient at that position
    Real  a_ij     = 0.0;
    Real  c_j      = 0.0;        ///< Objective coefficient

    // Saved values needed to reconstruct original solution
    Real  saved_lb  = 0.0;
    Real  saved_ub  = 0.0;
    Real  saved_rl  = 0.0;       ///< saved row_lower
    Real  saved_ru  = 0.0;       ///< saved row_upper
    Real  fix_value = 0.0;       ///< Value variable was fixed to
    Real  obj_offset = 0.0;      ///< Constant added to objective

    // For multi-variable operations (forcing row, probing)
    std::vector<Index> affected_cols;
    std::vector<Real>  saved_col_lbs;
    std::vector<Real>  saved_col_ubs;
};

/// LIFO undo stack.
class PresolveStack {
public:
    void push(PresolveRecord record) { records_.push_back(std::move(record)); }
    PresolveRecord pop() {
        auto r = std::move(records_.back());
        records_.pop_back();
        return r;
    }
    bool  empty() const { return records_.empty(); }
    Index size()  const { return static_cast<Index>(records_.size()); }

private:
    std::vector<PresolveRecord> records_;
};

// ── Presolve result ──────────────────────────────────────────────────────────

struct PresolveResult {
    Problem      reduced_problem;    ///< Smaller, tighter problem
    PresolveStack stack;             ///< Undo stack for postsolve
    bool         is_infeasible = false;
    bool         is_unbounded  = false;
    Index        rows_removed       = 0;
    Index        cols_removed       = 0;
    Index        bounds_tightened   = 0;
    Real         obj_offset         = 0.0;  ///< Constant added to objective during presolve
};

// ── Presolve engine ──────────────────────────────────────────────────────────

/// Applies problem reductions before solving.
/// Every reduction is reversible via the undo stack.
class Presolve {
public:
    /// Apply all enabled reduction rules until no more changes occur.
    /// Returns a PresolveResult with the reduced problem and undo stack.
    PresolveResult apply(const Problem& original);

    /// Recover full original-space solution from the reduced-space solution.
    Solution postsolve(const Solution& reduced_sol, const PresolveResult& result);

    // ── Configuration ─────────────────────────────────────────────────────
    void set_max_rounds(int n)          { max_rounds_ = n; }
    void set_enable_probing(bool b)     { enable_probing_ = b; }
    void set_log_level(LogLevel level)  { log_level_ = level; }

private:
    int      max_rounds_      = 20;
    bool     enable_probing_  = true;
    LogLevel log_level_       = LogLevel::INFO;

    // ── Reduction rules ────────────────────────────────────────────────────
    // Each returns true if the problem was changed.

    bool reduce_empty_rows         (Problem& p, PresolveStack& s, bool& infeasible);
    bool reduce_fixed_variables    (Problem& p, PresolveStack& s);
    bool reduce_singleton_rows     (Problem& p, PresolveStack& s, bool& infeasible);
    bool reduce_singleton_cols     (Problem& p, PresolveStack& s);
    bool reduce_forcing_rows       (Problem& p, PresolveStack& s, bool& infeasible);
    bool reduce_dominated_cols     (Problem& p, PresolveStack& s);
    bool reduce_duplicate_rows     (Problem& p, PresolveStack& s);
    bool reduce_duplicate_cols     (Problem& p, PresolveStack& s);
    bool tighten_implied_bounds    (Problem& p, PresolveStack& s, bool& infeasible);
    bool tighten_coefficients      (Problem& p, PresolveStack& s);  ///< MILP only
    bool run_probing               (Problem& p, PresolveStack& s, bool& infeasible);

    // ── Internal helpers ───────────────────────────────────────────────────
    void  compute_row_activities   (const Problem& p,
                                    std::vector<Real>& min_act,
                                    std::vector<Real>& max_act) const;
    bool  propagate_bounds         (Problem& p, bool& infeasible);

    /// Compact problem: renumber rows/cols after removals.
    Problem compact(const Problem& p,
                    const std::vector<bool>& row_active,
                    const std::vector<bool>& col_active,
                    std::vector<Index>& old_to_new_row,
                    std::vector<Index>& old_to_new_col) const;
};

} // namespace suplex
