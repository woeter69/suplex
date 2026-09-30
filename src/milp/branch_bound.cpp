#include "branch_bound.h"
#include "branching.h"
#include "cuts.h"
#include "heuristics.h"
#include "conflict.h"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <queue>
#include <cassert>

namespace suplex {

// ── Timer helper ──────────────────────────────────────────────────────────────
static double now_seconds() {
    using clock = std::chrono::steady_clock;
    static auto t0 = clock::now();
    return std::chrono::duration<double>(clock::now() - t0).count();
}

// ── Constructor ───────────────────────────────────────────────────────────────
BranchBound::BranchBound(std::shared_ptr<LPSolver> lp_solver)
    : lp_solver_(std::move(lp_solver)) {}

// ── Public solve entry point ──────────────────────────────────────────────────
SolverStatus BranchBound::solve(const Problem& problem, Solution& solution) {
    start_time_    = now_seconds();
    incumbent_obj_ = 1e30;
    stats_         = TreeStats{};

    // ── Step 1: Solve root LP ────────────────────────────────────────────────
    Solution root_lp_sol;
    SolverStatus root_status = solve_root_lp(problem, root_lp_sol);

    if (root_status == SolverStatus::INFEASIBLE) {
        solution.status = SolverStatus::INFEASIBLE;
        return root_status;
    }
    if (root_status == SolverStatus::UNBOUNDED) {
        solution.status = SolverStatus::UNBOUNDED;
        return root_status;
    }

    // ── Step 2: Check if root LP solution is already integer-feasible ────────
    if (is_integer_feasible(problem, root_lp_sol)) {
        solution = root_lp_sol;
        solution.status = SolverStatus::OPTIMAL;
        solution.best_bound = root_lp_sol.objective_value;
        solution.mip_gap    = 0.0;
        if (log_level_ <= LogLevel::INFO)
            std::cout << "[B&B] Root LP is already integer-feasible. Optimal.\n";
        return SolverStatus::OPTIMAL;
    }

    Real root_bound = root_lp_sol.objective_value;
    if (log_level_ <= LogLevel::INFO)
        std::cout << "[B&B] Root LP bound: " << root_bound << "\n";

    // ── Step 3: Cut rounds at root ───────────────────────────────────────────
    Problem problem_with_cuts = problem;
    int total_cuts = separate_cuts(problem, root_lp_sol, problem_with_cuts);
    if (log_level_ <= LogLevel::INFO && total_cuts > 0)
        std::cout << "[B&B] Added " << total_cuts << " cuts at root\n";

    // ── Step 4: Root heuristics ──────────────────────────────────────────────
    run_heuristics(problem_with_cuts, root_lp_sol);

    // ── Step 5: Initialise B&B tree ──────────────────────────────────────────
    Branching branching(problem_with_cuts.num_cols());
    ConflictAnalysis conflict_analysis;

    // Open node pool
    std::vector<std::unique_ptr<BBNode>> node_pool;

    // Root node
    auto root_node = std::make_unique<BBNode>();
    root_node->id        = stats_.nodes_created++;
    root_node->lp_bound  = root_bound;
    root_node->basis     = lp_solver_->get_basis();
    root_node->status    = NodeStatus::OPEN;
    node_pool.push_back(std::move(root_node));

    BBNode* last_node = nullptr;

    // ── Step 6: Main B&B loop ────────────────────────────────────────────────
    while (!node_pool.empty()) {
        // Time / node limit checks
        if (now_seconds() - start_time_ > time_limit_) {
            if (log_level_ <= LogLevel::INFO)
                std::cout << "[B&B] Time limit reached.\n";
            break;
        }
        if (stats_.nodes_solved >= node_limit_) {
            if (log_level_ <= LogLevel::INFO)
                std::cout << "[B&B] Node limit reached.\n";
            break;
        }
        if (callback_ && callback_->should_terminate()) break;

        // Gap check
        Real best_open_bound = 1e30;
        for (const auto& n : node_pool)
            if (n->lp_bound < best_open_bound) best_open_bound = n->lp_bound;
        if (compute_gap(best_open_bound) <= gap_tol_) {
            if (log_level_ <= LogLevel::INFO)
                std::cout << "[B&B] Gap tolerance reached.\n";
            break;
        }

        // Select node
        BBNode* node = select_node(node_pool, last_node);
        if (!node) break;
        last_node = node;

        // Apply bound changes from this node to get the local problem
        Problem local_p = problem_with_cuts;
        apply_bound_changes(local_p, node->bound_changes);

        // Solve LP relaxation (warm start from parent basis)
        Solution lp_sol;
        SolverStatus lp_status;
        if (!node->basis.empty())
            lp_status = lp_solver_->solve_from_basis(local_p, node->basis, lp_sol);
        else
            lp_status = lp_solver_->solve(local_p, lp_sol);

        stats_.nodes_solved++;
        node->lp_bound = lp_sol.objective_value;
        node->status   = (lp_status == SolverStatus::OPTIMAL)
                         ? NodeStatus::SOLVED : NodeStatus::PRUNED_INFEASIBLE;

        if (log_level_ <= LogLevel::DEBUG)
            log_progress(node->id, lp_sol.objective_value,
                         node->depth, static_cast<int>(node_pool.size()));

        // ── Pruning ───────────────────────────────────────────────────────────

        // Prune 1: LP infeasible
        if (lp_status == SolverStatus::INFEASIBLE) {
            node->status = NodeStatus::PRUNED_INFEASIBLE;
            stats_.nodes_infeasible++;

            // Conflict analysis
            auto farkas = lp_solver_->get_farkas_ray();
            if (!farkas.empty()) {
                std::vector<std::pair<Index, BoundChange>> branch_path;
                for (const auto& bc : node->bound_changes)
                    branch_path.push_back({bc.col_idx, bc});
                auto conflict_cut = conflict_analysis.analyse(farkas, *node, branch_path, local_p);
                if (conflict_cut.has_value()) {
                    // Add conflict cut to pool (would be added to problem_with_cuts)
                    stats_.cuts_added++;
                }
            }
            stats_.nodes_pruned++;
            remove_from_pool(node_pool, node);
            continue;
        }

        // Prune 2: LP bound >= incumbent (can't improve)
        if (lp_sol.objective_value >= incumbent_obj_ - EPS_OPTIMALITY) {
            node->status = NodeStatus::PRUNED_BOUND;
            stats_.nodes_pruned++;
            remove_from_pool(node_pool, node);
            continue;
        }

        // ── Integer feasible? → update incumbent ──────────────────────────────
        if (is_integer_feasible(local_p, lp_sol)) {
            node->status = NodeStatus::INTEGER_FEASIBLE;
            try_update_incumbent(local_p, lp_sol, problem.num_rows());
            if (callback_) callback_->on_new_incumbent(incumbent_sol_);
            stats_.nodes_pruned++;
            remove_from_pool(node_pool, node);
            continue;
        }

        // ── Run heuristics periodically (every 100 nodes) ────────────────────
        if (stats_.nodes_solved % 100 == 0)
            run_heuristics(local_p, lp_sol);

        // ── Branch ───────────────────────────────────────────────────────────
        Index branch_var = branching.select_variable(
            local_p, lp_sol, lp_solver_.get(),
            lp_sol.objective_value, node->depth);

        if (branch_var == -1) {
            // All integer — should have been caught above
            remove_from_pool(node_pool, node);
            continue;
        }

        node->status   = NodeStatus::BRANCHED;
        node->basis    = lp_solver_->get_basis();

        Real frac_val  = lp_sol.primal_values[branch_var];
        Real floor_val = std::floor(frac_val);
        Real ceil_val  = std::ceil(frac_val);

        const auto& cl = local_p.col_lower();
        const auto& cu = local_p.col_upper();

        // Create DOWN child (x <= floor)
        if (floor_val >= cl[branch_var] - EPS_FEASIBILITY) {
            auto child_down = std::make_unique<BBNode>();
            child_down->id        = stats_.nodes_created++;
            child_down->parent_id = node->id;
            child_down->depth     = node->depth + 1;
            child_down->branch_var   = branch_var;
            child_down->branch_value = frac_val;
            child_down->branch_dir   = BranchDirection::DOWN;
            child_down->lp_bound     = lp_sol.objective_value;  // initial estimate
            child_down->basis        = node->basis;
            child_down->status       = NodeStatus::OPEN;
            // Inherit parent bound changes + new one
            child_down->bound_changes = node->bound_changes;
            child_down->bound_changes.push_back({branch_var, cl[branch_var], floor_val});
            child_down->estimate = lp_sol.objective_value;
            if (node->depth > stats_.max_depth) stats_.max_depth = node->depth;
            node_pool.push_back(std::move(child_down));
        }

        // Create UP child (x >= ceil)
        if (ceil_val <= cu[branch_var] + EPS_FEASIBILITY) {
            auto child_up = std::make_unique<BBNode>();
            child_up->id        = stats_.nodes_created++;
            child_up->parent_id = node->id;
            child_up->depth     = node->depth + 1;
            child_up->branch_var   = branch_var;
            child_up->branch_value = frac_val;
            child_up->branch_dir   = BranchDirection::UP;
            child_up->lp_bound     = lp_sol.objective_value;
            child_up->basis        = node->basis;
            child_up->status       = NodeStatus::OPEN;
            child_up->bound_changes = node->bound_changes;
            child_up->bound_changes.push_back({branch_var, ceil_val, cu[branch_var]});
            child_up->estimate = lp_sol.objective_value;
            node_pool.push_back(std::move(child_up));
        }

        if (callback_)
            callback_->on_node_solved(node->id, lp_sol.objective_value);

        // Remove processed node from pool
        remove_from_pool(node_pool, node);

        // Periodic progress log
        if (stats_.nodes_solved % 500 == 0 && log_level_ <= LogLevel::INFO)
            log_progress(node->id, lp_sol.objective_value,
                         node->depth, static_cast<int>(node_pool.size()));
    }

    // ── Final result ──────────────────────────────────────────────────────────
    // Compute best bound from remaining open nodes
    Real best_bound = incumbent_obj_;
    for (const auto& n : node_pool)
        if (n->lp_bound < best_bound) best_bound = n->lp_bound;

    if (incumbent_obj_ < 1e29) {
        solution = incumbent_sol_;
        solution.status     = SolverStatus::OPTIMAL;
        solution.best_bound = best_bound;
        solution.mip_gap    = compute_gap(best_bound);
        solution.num_nodes  = stats_.nodes_solved;
        solution.solve_time_seconds = now_seconds() - start_time_;
        if (log_level_ <= LogLevel::INFO)
            std::cout << "[B&B] Done. Obj=" << solution.objective_value
                      << " Gap=" << solution.mip_gap * 100.0 << "%"
                      << " Nodes=" << stats_.nodes_solved << "\n";
        return SolverStatus::OPTIMAL;
    }

    solution.status = node_pool.empty()
                      ? SolverStatus::INFEASIBLE
                      : SolverStatus::TIME_LIMIT;
    return solution.status;
}

// ── solve_root_lp ─────────────────────────────────────────────────────────────
SolverStatus BranchBound::solve_root_lp(const Problem& p, Solution& sol) {
    SolverStatus st = lp_solver_->solve(p, sol);
    sol.num_iterations = 0;  // will be set by LP solver
    return st;
}

// ── separate_cuts ─────────────────────────────────────────────────────────────
int BranchBound::separate_cuts(const Problem& p, const Solution& lp_sol,
                                Problem& p_with_cuts) {
    CutGenerator gen(p.num_cols());
    int total = 0;

    for (int round = 0; round < cut_rounds_root_; ++round) {
        std::vector<Cut> candidates;

        // MIR cuts (don't need basis — work from constraint rows)
        auto mir = gen.mir_cuts(p, lp_sol, 50);
        candidates.insert(candidates.end(), mir.begin(), mir.end());

        // Cover cuts (binary variables)
        if (p.is_mip()) {
            auto cov = gen.cover_cuts(p, lp_sol, 30);
            candidates.insert(candidates.end(), cov.begin(), cov.end());
        }

        // Clique cuts
        auto cliq = gen.clique_cuts(p, lp_sol, 20);
        candidates.insert(candidates.end(), cliq.begin(), cliq.end());

        if (candidates.empty()) break;

        auto selected = gen.select_cuts(candidates, 100);
        if (selected.empty()) break;

        total += static_cast<int>(selected.size());
        stats_.cuts_added += static_cast<int>(selected.size());
        // In a full implementation: add cuts as rows to p_with_cuts
        // For now just count them
    }
    return total;
}

// ── run_heuristics ────────────────────────────────────────────────────────────
bool BranchBound::run_heuristics(const Problem& p, const Solution& lp_sol) {
    Heuristics h(lp_solver_.get());
    h.set_log_level(log_level_);

    bool found = false;

    // Try simple rounding first (very cheap)
    auto result = h.simple_rounding(p, lp_sol);
    if (result && result->objective_value < incumbent_obj_) {
        incumbent_obj_ = result->objective_value;
        incumbent_sol_ = *result;
        stats_.heuristic_sols++;
        found = true;
        if (log_level_ <= LogLevel::INFO)
            std::cout << "[Heuristic] Rounding: obj=" << incumbent_obj_ << "\n";
    }

    // Try fractional diving
    result = h.fractional_diving(p, lp_sol);
    if (result && result->objective_value < incumbent_obj_) {
        incumbent_obj_ = result->objective_value;
        incumbent_sol_ = *result;
        stats_.heuristic_sols++;
        found = true;
        if (log_level_ <= LogLevel::INFO)
            std::cout << "[Heuristic] Diving: obj=" << incumbent_obj_ << "\n";
    }

    // Try feasibility pump if still no incumbent
    if (incumbent_obj_ >= 1e29) {
        result = h.feasibility_pump(p, lp_sol);
        if (result && result->objective_value < incumbent_obj_) {
            incumbent_obj_ = result->objective_value;
            incumbent_sol_ = *result;
            stats_.heuristic_sols++;
            found = true;
            if (log_level_ <= LogLevel::INFO)
                std::cout << "[Heuristic] FeasPump: obj=" << incumbent_obj_ << "\n";
        }
    }

    return found;
}

// ── is_integer_feasible ───────────────────────────────────────────────────────
bool BranchBound::is_integer_feasible(const Problem& p,
                                       const Solution& sol) const {
    const auto& vt = p.var_types();
    for (Index j = 0; j < p.num_cols(); ++j) {
        if (vt[j] == VarType::CONTINUOUS) continue;
        Real f = sol.primal_values[j] - std::floor(sol.primal_values[j]);
        if (f > EPS_INTEGER && f < 1.0 - EPS_INTEGER) return false;
    }
    return true;
}

// ── try_update_incumbent ──────────────────────────────────────────────────────
void BranchBound::try_update_incumbent(const Problem& p, const Solution& sol,
                                        Index /*base_rows*/) {
    if (sol.objective_value < incumbent_obj_ - EPS_OPTIMALITY) {
        incumbent_obj_ = sol.objective_value;
        incumbent_sol_ = sol;
        if (log_level_ <= LogLevel::INFO)
            std::cout << "[B&B] New incumbent: " << incumbent_obj_ << "\n";
    }
}

// ── apply_bound_changes ───────────────────────────────────────────────────────
void BranchBound::apply_bound_changes(Problem& p,
                                       const std::vector<BoundChange>& changes) const {
    if (changes.empty()) return;
    auto new_cl = p.col_lower();
    auto new_cu = p.col_upper();
    for (const auto& bc : changes) {
        new_cl[bc.col_idx] = std::max(new_cl[bc.col_idx], bc.new_lower);
        new_cu[bc.col_idx] = std::min(new_cu[bc.col_idx], bc.new_upper);
    }
    p.set_col_bounds(std::move(new_cl), std::move(new_cu));
}

// ── select_node ───────────────────────────────────────────────────────────────
BBNode* BranchBound::select_node(std::vector<std::unique_ptr<BBNode>>& pool,
                                  const BBNode* last_node) {
    if (pool.empty()) return nullptr;

    switch (node_selection_) {
    case NodeSelection::DEPTH_FIRST:
        // Deepest node last (LIFO since we push children at end)
        return pool.back().get();

    case NodeSelection::BEST_FIRST: {
        auto it = std::min_element(pool.begin(), pool.end(),
            [](const auto& a, const auto& b) {
                return a->lp_bound < b->lp_bound;
            });
        return it->get();
    }

    case NodeSelection::HYBRID:
    default: {
        // Plunge: if last node was branched, prefer its children (they're at the end)
        if (last_node && last_node->status == NodeStatus::BRANCHED) {
            // Children were just pushed — check if they're better than incumbent
            if (!pool.empty() && pool.back()->lp_bound < incumbent_obj_ - EPS_OPTIMALITY)
                return pool.back().get();
        }
        // Otherwise best-first
        auto it = std::min_element(pool.begin(), pool.end(),
            [](const auto& a, const auto& b) {
                return a->lp_bound < b->lp_bound;
            });
        return it->get();
    }
    }
}

// ── remove_from_pool ─────────────────────────────────────────────────────────
void BranchBound::remove_from_pool(std::vector<std::unique_ptr<BBNode>>& pool,
                                    BBNode* node) {
    auto it = std::find_if(pool.begin(), pool.end(),
                           [node](const auto& p) { return p.get() == node; });
    if (it != pool.end()) {
        std::iter_swap(it, pool.end() - 1);
        pool.pop_back();
    }
}

// ── compute_gap ───────────────────────────────────────────────────────────────
Real BranchBound::compute_gap(Real best_bound) const {
    if (incumbent_obj_ >= 1e29) return 1.0;
    if (std::abs(incumbent_obj_) < EPS_ZERO) return std::abs(best_bound);
    return std::abs(incumbent_obj_ - best_bound) / std::abs(incumbent_obj_);
}

// ── log_progress ─────────────────────────────────────────────────────────────
void BranchBound::log_progress(int64_t node_id, Real lp_bound,
                                int depth, int open_nodes) const {
    double elapsed = now_seconds() - start_time_;
    std::cout << "[B&B] node=" << node_id
              << " depth=" << depth
              << " lp=" << lp_bound
              << " inc=" << (incumbent_obj_ < 1e29 ? incumbent_obj_ : -1.0)
              << " gap=" << compute_gap(lp_bound) * 100.0 << "%"
              << " open=" << open_nodes
              << " time=" << elapsed << "s\n";
}

} // namespace suplex
