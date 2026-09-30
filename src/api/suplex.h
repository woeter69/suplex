#ifndef SUPLEX_C_API_H
#define SUPLEX_C_API_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SuplexSolver SuplexSolver;

enum SuplexStatus {
    SUPLEX_STATUS_NOT_STARTED = 0,
    SUPLEX_STATUS_OPTIMAL = 1,
    SUPLEX_STATUS_INFEASIBLE = 2,
    SUPLEX_STATUS_UNBOUNDED = 3,
    SUPLEX_STATUS_INF_OR_UNBD = 4,
    SUPLEX_STATUS_ITERATION_LIMIT = 5,
    SUPLEX_STATUS_TIME_LIMIT = 6,
    SUPLEX_STATUS_NUMERICAL_ERROR = 7,
    SUPLEX_STATUS_USER_INTERRUPT = 8
};

enum SuplexAlgorithm {
    SUPLEX_ALGORITHM_AUTO = 0,
    SUPLEX_ALGORITHM_PRIMAL = 1,
    SUPLEX_ALGORITHM_DUAL = 2,
    SUPLEX_ALGORITHM_IPM = 3
};

enum SuplexVariableType {
    SUPLEX_VARIABLE_CONTINUOUS = 0,
    SUPLEX_VARIABLE_INTEGER = 1,
    SUPLEX_VARIABLE_BINARY = 2
};

enum SuplexErrorCode {
    SUPLEX_OK = 0,
    SUPLEX_ERROR_INVALID_ARGUMENT = -1,
    SUPLEX_ERROR_IO = -2,
    SUPLEX_ERROR_INVALID_PROBLEM = -3,
    SUPLEX_ERROR_BUFFER_TOO_SMALL = -4
};

/** Allocate a solver handle, or return NULL if allocation fails. */
SuplexSolver* suplex_create(void);
/** Release a handle. Passing NULL is valid. */
void suplex_destroy(SuplexSolver* solver);
/** Load an MPS model. Returns SUPLEX_OK or a negative SuplexErrorCode. */
int suplex_read_mps(SuplexSolver* solver, const char* filename);
/** Load a CPLEX LP model. Returns SUPLEX_OK or a negative SuplexErrorCode. */
int suplex_read_lp(SuplexSolver* solver, const char* filename);
/** Reset a programmatic model to the requested dimensions. */
int suplex_set_dimensions(SuplexSolver* solver, int rows, int cols);
/** Set dense objective coefficients and 0=minimize/1=maximize sense. */
int suplex_set_objective(SuplexSolver* solver, const double* coeffs, int sense);
/** Set a zero-based coordinate-format constraint matrix. */
int suplex_set_constraint_matrix(SuplexSolver* solver, int nnz, const int* row_indices,
                                 const int* col_indices, const double* values);
/** Set dense row lower and upper bounds. */
int suplex_set_row_bounds(SuplexSolver* solver, const double* lower, const double* upper);
/** Set dense column lower and upper bounds. */
int suplex_set_col_bounds(SuplexSolver* solver, const double* lower, const double* upper);
/** Set variable types using SuplexVariableType values. */
int suplex_set_var_types(SuplexSolver* solver, const int* types);
/** Select a SuplexAlgorithm. */
int suplex_set_algorithm(SuplexSolver* solver, int algorithm);
int suplex_set_presolve(SuplexSolver* solver, int enabled);
int suplex_set_time_limit(SuplexSolver* solver, double seconds);
int suplex_set_log_level(SuplexSolver* solver, int level);
int suplex_set_threads(SuplexSolver* solver, int threads);
int suplex_set_gpu(SuplexSolver* solver, int enabled);
int suplex_set_mip_gap(SuplexSolver* solver, double gap);
/** Solve and return a nonnegative SuplexStatus, or a negative error code. */
int suplex_solve(SuplexSolver* solver);
/** Return the current SuplexStatus, or -1 for a NULL handle. */
int suplex_get_status(const SuplexSolver* solver);
double suplex_get_objective(const SuplexSolver* solver);
/** Copy primal values; returns copied count or a negative error code. */
int suplex_get_solution(const SuplexSolver* solver, double* values, int size);
int suplex_get_dual_values(const SuplexSolver* solver, double* values, int size);
int suplex_get_reduced_costs(const SuplexSolver* solver, double* values, int size);
double suplex_get_solve_time(const SuplexSolver* solver);
int64_t suplex_get_iterations(const SuplexSolver* solver);
int64_t suplex_get_nodes(const SuplexSolver* solver);
double suplex_get_mip_gap(const SuplexSolver* solver);
const char* suplex_get_error(const SuplexSolver* solver);
const char* suplex_version(void);

#ifdef __cplusplus
}
#endif
#endif
