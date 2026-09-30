#include "suplex_solver.h"

#include "src/io/lp_parser.h"
#include "src/io/mps_parser.h"
#include "src/milp/milp_solver.h"
#include "src/simplex/dual_simplex.h"
#include "src/simplex/primal_simplex.h"

#include <chrono>
#include <memory>
#include <utility>

namespace suplex {

bool Suplex::read_mps(const std::string& filename) {
    MPSReader reader;
    if (!reader.read(filename, problem_)) { error_ = reader.last_error(); return false; }
    column_names_ = reader.column_names();
    row_names_ = reader.row_names();
    error_.clear();
    return true;
}

bool Suplex::read_lp(const std::string& filename) {
    LPReader reader;
    if (!reader.read(filename, problem_)) { error_ = reader.last_error(); return false; }
    column_names_ = reader.column_names();
    row_names_ = reader.row_names();
    error_.clear();
    return true;
}

void Suplex::set_problem(Problem problem) {
    problem_ = std::move(problem);
    column_names_.clear();
    row_names_.clear();
    error_.clear();
}

SolverStatus Suplex::solve() {
    solution_.reset();
    error_.clear();
    if (!problem_.validate()) {
        error_ = "problem data is dimensionally invalid or has inconsistent bounds";
        solution_.status = SolverStatus::NUMERICAL_ERROR;
        return solution_.status;
    }
    if (problem_.is_qp()) {
        error_ = "QP solving is unavailable because the P2 IPM/QP module is not present";
        solution_.status = SolverStatus::NUMERICAL_ERROR;
        return solution_.status;
    }
    if (algorithm_ == Algorithm::IPM) {
        error_ = "IPM was requested but the P2 IPM module is not present";
        solution_.status = SolverStatus::NUMERICAL_ERROR;
        return solution_.status;
    }
    (void)gpu_;
    (void)threads_;
    const auto start = std::chrono::steady_clock::now();

    std::shared_ptr<LPSolver> lp_solver;
    if (algorithm_ == Algorithm::PRIMAL_SIMPLEX) lp_solver = std::make_shared<PrimalSimplexSolver>();
    else lp_solver = std::make_shared<DualSimplexSolver>();
    lp_solver->set_time_limit(time_limit_);
    lp_solver->set_log_level(log_level_);

    SolverStatus status;
    if (problem_.is_mip()) {
        MILPSolver solver;
        solver.set_lp_solver(lp_solver);
        // The inherited postsolve implementation cannot currently reconstruct
        // eliminated columns reliably, including the all-binary singleton case.
        // Keep the public option ABI-stable but preserve correctness by solving
        // in original space until P3 completes that mapping.
        (void)presolve_;
        solver.set_presolve(false);
        solver.set_gap_tolerance(mip_gap_);
        solver.set_time_limit(time_limit_);
        solver.set_log_level(log_level_);
        status = solver.solve(problem_, solution_);
    } else {
        // P3's current postsolve stack does not retain a complete original-to-
        // reduced column map. Until that contract is completed, solving an LP
        // in original space is required to preserve objective/primal values.
        // MILP is likewise kept in original space above.
        (void)presolve_;
        status = lp_solver->solve(problem_, solution_);
    }

    solution_.solve_time_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    return status;
}

} // namespace suplex
