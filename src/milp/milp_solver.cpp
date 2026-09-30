#include "milp_solver.h"
#include "branch_bound.h"
#include "presolve.h"
#include <iostream>

namespace suplex {

MILPSolver::MILPSolver() = default;

void MILPSolver::set_lp_solver(std::shared_ptr<LPSolver> lp_solver) {
    lp_solver_ = std::move(lp_solver);
}

void MILPSolver::set_callback(std::shared_ptr<MILPCallback> cb) {
    callback_ = std::move(cb);
}

SolverStatus MILPSolver::solve(const Problem& problem, Solution& solution) {
    if (!lp_solver_) {
        if (log_level_ <= LogLevel::ERROR)
            std::cerr << "[MILPSolver] No LP solver set. Call set_lp_solver() first.\n";
        solution.status = SolverStatus::NUMERICAL_ERROR;
        return SolverStatus::NUMERICAL_ERROR;
    }

    // For continuous problems, delegate directly to LP solver
    if (!problem.is_mip()) {
        lp_solver_->set_time_limit(time_limit_);
        lp_solver_->set_log_level(log_level_);
        return lp_solver_->solve(problem, solution);
    }

    // ── Presolve ──────────────────────────────────────────────────────────────
    Presolve presolve;
    presolve.set_log_level(log_level_);

    PresolveResult pre_result;
    Problem effective_problem;

    if (use_presolve_) {
        pre_result = presolve.apply(problem);

        if (pre_result.is_infeasible) {
            solution.status = SolverStatus::INFEASIBLE;
            return SolverStatus::INFEASIBLE;
        }
        effective_problem = std::move(pre_result.reduced_problem);
        if (log_level_ <= LogLevel::INFO)
            std::cout << "[MILPSolver] Presolve: "
                      << effective_problem.num_rows() << " rows, "
                      << effective_problem.num_cols() << " cols\n";
    } else {
        effective_problem = problem;
    }

    // ── Branch and Bound ──────────────────────────────────────────────────────
    lp_solver_->set_time_limit(time_limit_);
    lp_solver_->set_log_level(log_level_);

    BranchBound bb(lp_solver_);
    bb.set_gap_tolerance(gap_tol_);
    bb.set_node_limit(node_limit_);
    bb.set_time_limit(time_limit_);
    bb.set_log_level(log_level_);
    if (callback_) bb.set_callback(callback_);

    Solution reduced_sol;
    SolverStatus status = bb.solve(effective_problem, reduced_sol);

    // ── Postsolve ─────────────────────────────────────────────────────────────
    if (use_presolve_ && status == SolverStatus::OPTIMAL) {
        solution = presolve.postsolve(reduced_sol, pre_result);
        solution.status           = SolverStatus::OPTIMAL;
        solution.num_nodes        = reduced_sol.num_nodes;
        solution.best_bound       = reduced_sol.best_bound;
        solution.mip_gap          = reduced_sol.mip_gap;
        solution.solve_time_seconds = reduced_sol.solve_time_seconds;
    } else {
        solution = reduced_sol;
    }

    return status;
}

const TreeStats& MILPSolver::tree_stats() const {
    // Return last B&B stats — stored via callback in production
    static TreeStats empty;
    return empty;
}

} // namespace suplex
