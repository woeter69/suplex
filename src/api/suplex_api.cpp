#include "suplex.h"

#include "suplex_solver.h"
#include "src/core/sparse_matrix.h"

#include <algorithm>
#include <new>
#include <string>
#include <vector>

struct SuplexSolver {
    suplex::Suplex facade;
    suplex::Problem problem;
    std::string error;
    bool programmatic = false;
};

namespace {

int invalid(SuplexSolver* solver, const char* message) {
    if (solver) solver->error = message;
    return SUPLEX_ERROR_INVALID_ARGUMENT;
}

template <typename T>
int copy_result(const SuplexSolver* solver, const std::vector<T>& source, double* values, int size) {
    if (!solver || !values || size < 0) return SUPLEX_ERROR_INVALID_ARGUMENT;
    if (size < static_cast<int>(source.size())) return SUPLEX_ERROR_BUFFER_TOO_SMALL;
    std::copy(source.begin(), source.end(), values);
    return static_cast<int>(source.size());
}

} // namespace

extern "C" {

SuplexSolver* suplex_create(void) {
    return new (std::nothrow) SuplexSolver;
}

void suplex_destroy(SuplexSolver* solver) { delete solver; }

int suplex_read_mps(SuplexSolver* solver, const char* filename) {
    if (!solver || !filename) return invalid(solver, "solver and filename are required");
    if (!solver->facade.read_mps(filename)) { solver->error = solver->facade.last_error(); return SUPLEX_ERROR_IO; }
    solver->programmatic = false; solver->error.clear(); return SUPLEX_OK;
}

int suplex_read_lp(SuplexSolver* solver, const char* filename) {
    if (!solver || !filename) return invalid(solver, "solver and filename are required");
    if (!solver->facade.read_lp(filename)) { solver->error = solver->facade.last_error(); return SUPLEX_ERROR_IO; }
    solver->programmatic = false; solver->error.clear(); return SUPLEX_OK;
}

int suplex_set_dimensions(SuplexSolver* solver, int rows, int cols) {
    if (!solver || rows < 0 || cols < 0) return invalid(solver, "dimensions must be nonnegative");
    solver->problem.set_dimensions(rows, cols); solver->programmatic = true; solver->error.clear(); return SUPLEX_OK;
}

int suplex_set_objective(SuplexSolver* solver, const double* coeffs, int sense) {
    if (!solver || (!coeffs && solver->problem.num_cols() != 0) || (sense != 0 && sense != 1))
        return invalid(solver, "invalid objective data or sense");
    std::vector<double> objective;
    if (solver->problem.num_cols() > 0) objective.assign(coeffs, coeffs + solver->problem.num_cols());
    solver->problem.set_objective(std::move(objective));
    solver->problem.set_objective_sense(sense == 0 ? suplex::ObjectiveSense::MINIMIZE : suplex::ObjectiveSense::MAXIMIZE);
    solver->programmatic = true; return SUPLEX_OK;
}

int suplex_set_constraint_matrix(SuplexSolver* solver, int nnz, const int* rows,
                                 const int* cols, const double* values) {
    if (!solver || nnz < 0 || (nnz > 0 && (!rows || !cols || !values))) return invalid(solver, "invalid matrix arrays");
    std::vector<suplex::Index> r, c;
    std::vector<double> v;
    r.reserve(static_cast<std::size_t>(nnz)); c.reserve(static_cast<std::size_t>(nnz)); v.reserve(static_cast<std::size_t>(nnz));
    for (int k = 0; k < nnz; ++k) {
        if (rows[k] < 0 || rows[k] >= solver->problem.num_rows() || cols[k] < 0 || cols[k] >= solver->problem.num_cols())
            return invalid(solver, "matrix index is out of range");
        r.push_back(rows[k]); c.push_back(cols[k]); v.push_back(values[k]);
    }
    solver->problem.set_constraint_matrix(suplex::SparseMatrixCSC::from_triplets(
        solver->problem.num_rows(), solver->problem.num_cols(), r, c, v));
    solver->programmatic = true; return SUPLEX_OK;
}

int suplex_set_row_bounds(SuplexSolver* solver, const double* lower, const double* upper) {
    if (!solver || ((!lower || !upper) && solver->problem.num_rows() != 0)) return invalid(solver, "row bound arrays are required");
    std::vector<double> lower_values, upper_values;
    if (solver->problem.num_rows() > 0) {
        lower_values.assign(lower, lower + solver->problem.num_rows());
        upper_values.assign(upper, upper + solver->problem.num_rows());
    }
    solver->problem.set_row_bounds(std::move(lower_values), std::move(upper_values));
    solver->programmatic = true; return SUPLEX_OK;
}

int suplex_set_col_bounds(SuplexSolver* solver, const double* lower, const double* upper) {
    if (!solver || ((!lower || !upper) && solver->problem.num_cols() != 0)) return invalid(solver, "column bound arrays are required");
    std::vector<double> lower_values, upper_values;
    if (solver->problem.num_cols() > 0) {
        lower_values.assign(lower, lower + solver->problem.num_cols());
        upper_values.assign(upper, upper + solver->problem.num_cols());
    }
    solver->problem.set_col_bounds(std::move(lower_values), std::move(upper_values));
    solver->programmatic = true; return SUPLEX_OK;
}

int suplex_set_var_types(SuplexSolver* solver, const int* types) {
    if (!solver || (!types && solver->problem.num_cols() != 0)) return invalid(solver, "variable type array is required");
    std::vector<suplex::VarType> result;
    result.reserve(static_cast<std::size_t>(solver->problem.num_cols()));
    for (suplex::Index j = 0; j < solver->problem.num_cols(); ++j) {
        if (types[j] < 0 || types[j] > 2) return invalid(solver, "variable types must be 0, 1, or 2");
        result.push_back(static_cast<suplex::VarType>(types[j]));
    }
    solver->problem.set_var_types(std::move(result)); solver->programmatic = true; return SUPLEX_OK;
}

int suplex_set_algorithm(SuplexSolver* solver, int value) {
    if (!solver || value < 0 || value > 3) return invalid(solver, "algorithm must be in [0,3]");
    solver->facade.set_algorithm(static_cast<suplex::Algorithm>(value)); return SUPLEX_OK;
}
int suplex_set_presolve(SuplexSolver* solver, int value) { if (!solver) return SUPLEX_ERROR_INVALID_ARGUMENT; solver->facade.set_presolve(value != 0); return SUPLEX_OK; }
int suplex_set_time_limit(SuplexSolver* solver, double value) { if (!solver || value <= 0.0) return invalid(solver, "time limit must be positive"); solver->facade.set_time_limit(value); return SUPLEX_OK; }
int suplex_set_log_level(SuplexSolver* solver, int value) { if (!solver || value < 0 || value > 5) return invalid(solver, "log level must be in [0,5]"); solver->facade.set_log_level(static_cast<suplex::LogLevel>(value)); return SUPLEX_OK; }
int suplex_set_threads(SuplexSolver* solver, int value) { if (!solver || value < 0) return invalid(solver, "thread count must be nonnegative"); solver->facade.set_threads(value); return SUPLEX_OK; }
int suplex_set_gpu(SuplexSolver* solver, int value) { if (!solver) return SUPLEX_ERROR_INVALID_ARGUMENT; solver->facade.set_gpu(value != 0); return SUPLEX_OK; }
int suplex_set_mip_gap(SuplexSolver* solver, double value) { if (!solver || value < 0.0) return invalid(solver, "MIP gap must be nonnegative"); solver->facade.set_mip_gap(value); return SUPLEX_OK; }

int suplex_solve(SuplexSolver* solver) {
    if (!solver) return SUPLEX_ERROR_INVALID_ARGUMENT;
    if (solver->programmatic) solver->facade.set_problem(solver->problem);
    const auto status = solver->facade.solve();
    solver->error = solver->facade.last_error();
    return static_cast<int>(status);
}

int suplex_get_status(const SuplexSolver* solver) { return solver ? static_cast<int>(solver->facade.solution().status) : -1; }
double suplex_get_objective(const SuplexSolver* solver) { return solver ? solver->facade.solution().objective_value : 0.0; }
int suplex_get_solution(const SuplexSolver* solver, double* values, int size) { return solver ? copy_result(solver, solver->facade.solution().primal_values, values, size) : SUPLEX_ERROR_INVALID_ARGUMENT; }
int suplex_get_dual_values(const SuplexSolver* solver, double* values, int size) { return solver ? copy_result(solver, solver->facade.solution().dual_values, values, size) : SUPLEX_ERROR_INVALID_ARGUMENT; }
int suplex_get_reduced_costs(const SuplexSolver* solver, double* values, int size) { return solver ? copy_result(solver, solver->facade.solution().reduced_costs, values, size) : SUPLEX_ERROR_INVALID_ARGUMENT; }
double suplex_get_solve_time(const SuplexSolver* solver) { return solver ? solver->facade.solution().solve_time_seconds : 0.0; }
int64_t suplex_get_iterations(const SuplexSolver* solver) { return solver ? solver->facade.solution().num_iterations : 0; }
int64_t suplex_get_nodes(const SuplexSolver* solver) { return solver ? solver->facade.solution().num_nodes : 0; }
double suplex_get_mip_gap(const SuplexSolver* solver) { return solver ? solver->facade.solution().mip_gap : 0.0; }
const char* suplex_get_error(const SuplexSolver* solver) { return solver ? solver->error.c_str() : "null solver"; }
const char* suplex_version(void) { return "0.1.0"; }

} // extern "C"
