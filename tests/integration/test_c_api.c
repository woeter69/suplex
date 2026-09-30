#include "src/api/suplex.h"

#include <math.h>
#include <stddef.h>

int main(void) {
    SuplexSolver* solver = suplex_create();
    if (solver == NULL) return 1;
    const double objective[] = {1.0};
    const int rows[] = {0};
    const int columns[] = {0};
    const double coefficients[] = {1.0};
    const double row_lower[] = {2.0};
    const double row_upper[] = {1e30};
    const double col_lower[] = {0.0};
    const double col_upper[] = {10.0};
    if (suplex_set_dimensions(solver, 1, 1) != SUPLEX_OK ||
        suplex_set_objective(solver, objective, 0) != SUPLEX_OK ||
        suplex_set_constraint_matrix(solver, 1, rows, columns, coefficients) != SUPLEX_OK ||
        suplex_set_row_bounds(solver, row_lower, row_upper) != SUPLEX_OK ||
        suplex_set_col_bounds(solver, col_lower, col_upper) != SUPLEX_OK ||
        suplex_set_algorithm(solver, SUPLEX_ALGORITHM_PRIMAL) != SUPLEX_OK ||
        suplex_set_presolve(solver, 0) != SUPLEX_OK ||
        suplex_solve(solver) != SUPLEX_STATUS_OPTIMAL) {
        suplex_destroy(solver);
        return 2;
    }
    double value = 0.0;
    const int copied = suplex_get_solution(solver, &value, 1);
    suplex_destroy(solver);
    return copied == 1 && fabs(value - 2.0) <= 1e-8 ? 0 : 3;
}
