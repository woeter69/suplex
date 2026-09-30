#include "tests/unit/test_framework.h"

#include "src/api/suplex.h"
#include "src/api/suplex_solver.h"
#include "src/gpu/gpu_manager.h"
#include "src/io/lp_parser.h"
#include "src/io/mps_parser.h"
#include "src/io/sol_writer.h"

#include <filesystem>
#include <vector>

using namespace suplex;

namespace {
std::string fixture(const char* name) {
    return (std::filesystem::path(SUPLEX_SOURCE_DIR) / "data" / "examples" / name).string();
}

void test_mps_and_lp_equivalence() {
    MPSReader mps; LPReader lp; Problem a, b;
    TEST_ASSERT(mps.read(fixture("tiny.mps"), a));
    TEST_ASSERT(lp.read(fixture("tiny.lp"), b));
    TEST_ASSERT(a.validate() && b.validate());
    TEST_ASSERT(a.num_rows() == b.num_rows());
    TEST_ASSERT(a.num_cols() == b.num_cols());
    TEST_ASSERT(a.constraint_matrix().nnz == b.constraint_matrix().nnz);
    TEST_ASSERT(mps.column_names().size() == 2);
    TEST_ASSERT(lp.row_names().size() == 2);
}
REGISTER_TEST(test_mps_and_lp_equivalence);

void test_mps_features() {
    MPSReader reader; Problem problem;
    TEST_ASSERT(reader.read(fixture("features.mps"), problem));
    TEST_ASSERT(problem.sense() == ObjectiveSense::MAXIMIZE);
    TEST_ASSERT(problem.num_cols() == 2 && problem.num_rows() == 2);
    TEST_ASSERT(problem.var_types()[0] == VarType::CONTINUOUS);
    TEST_ASSERT(problem.var_types()[1] == VarType::BINARY);
    TEST_ASSERT(problem.col_lower()[0] == -INF && problem.col_upper()[0] == INF);
    TEST_ASSERT(problem.col_lower()[1] == 0.0 && problem.col_upper()[1] == 1.0);
    TEST_ASSERT_NEAR(problem.row_lower()[0], 6.0, 1e-12);
    TEST_ASSERT_NEAR(problem.row_upper()[0], 10.0, 1e-12);
    TEST_ASSERT_NEAR(problem.row_lower()[1], 1.0, 1e-12);
    TEST_ASSERT_NEAR(problem.row_upper()[1], 3.0, 1e-12);
}
REGISTER_TEST(test_mps_features);

void test_lp_compact_bounds_and_types() {
    LPReader reader; Problem problem;
    TEST_ASSERT(reader.read(fixture("compact_bounds.lp"), problem));
    TEST_ASSERT(problem.sense() == ObjectiveSense::MAXIMIZE);
    TEST_ASSERT(problem.num_cols() == 2);
    TEST_ASSERT_NEAR(problem.objective()[0], 2.0, 1e-12);
    TEST_ASSERT_NEAR(problem.objective()[1], -1.0, 1e-12);
    TEST_ASSERT(problem.var_types()[0] == VarType::INTEGER);
    TEST_ASSERT(problem.col_lower()[0] == 0.0 && problem.col_upper()[0] == 4.0);
    TEST_ASSERT(problem.col_lower()[1] == -INF && problem.col_upper()[1] == INF);
}
REGISTER_TEST(test_lp_compact_bounds_and_types);

void test_parser_error_diagnostic() {
    MPSReader reader; Problem problem;
    TEST_ASSERT(!reader.read(fixture("invalid_missing_end.mps"), problem));
    TEST_ASSERT(reader.last_error().find("missing ENDATA") != std::string::npos);
    TEST_ASSERT(reader.last_error().find("line") != std::string::npos);
}
REGISTER_TEST(test_parser_error_diagnostic);

void test_top_level_file_solve() {
    Suplex solver;
    TEST_ASSERT(solver.read_mps(fixture("tiny.mps")));
    solver.set_algorithm(Algorithm::PRIMAL_SIMPLEX);
    solver.set_presolve(false);
    TEST_ASSERT(solver.solve() == SolverStatus::OPTIMAL);
    TEST_ASSERT_NEAR(solver.solution().objective_value, 3.0, 1e-7);
    TEST_ASSERT(solver.solution().primal_values.size() == 2);
}
REGISTER_TEST(test_top_level_file_solve);

void test_c_api_programmatic_problem() {
    SuplexSolver* solver = suplex_create();
    TEST_ASSERT(solver != nullptr);
    TEST_ASSERT(suplex_set_dimensions(solver, 1, 1) == SUPLEX_OK);
    const double objective[] = {1.0};
    const int row[] = {0}, col[] = {0};
    const double coefficient[] = {1.0};
    const double row_lower[] = {2.0}, row_upper[] = {INF};
    const double col_lower[] = {0.0}, col_upper[] = {10.0};
    TEST_ASSERT(suplex_set_objective(solver, objective, 0) == SUPLEX_OK);
    TEST_ASSERT(suplex_set_constraint_matrix(solver, 1, row, col, coefficient) == SUPLEX_OK);
    TEST_ASSERT(suplex_set_row_bounds(solver, row_lower, row_upper) == SUPLEX_OK);
    TEST_ASSERT(suplex_set_col_bounds(solver, col_lower, col_upper) == SUPLEX_OK);
    TEST_ASSERT(suplex_set_algorithm(solver, 1) == SUPLEX_OK);
    TEST_ASSERT(suplex_set_presolve(solver, 0) == SUPLEX_OK);
    TEST_ASSERT(suplex_solve(solver) == static_cast<int>(SolverStatus::OPTIMAL));
    double value = 0.0;
    TEST_ASSERT(suplex_get_solution(solver, &value, 1) == 1);
    TEST_ASSERT_NEAR(value, 2.0, 1e-7);
    suplex_destroy(solver);
}
REGISTER_TEST(test_c_api_programmatic_problem);

void test_c_api_validation_and_buffer_size() {
    SuplexSolver* solver = suplex_create();
    TEST_ASSERT(solver != nullptr);
    TEST_ASSERT(suplex_set_dimensions(solver, -1, 2) == SUPLEX_ERROR_INVALID_ARGUMENT);
    TEST_ASSERT(std::string(suplex_get_error(solver)).find("nonnegative") != std::string::npos);
    TEST_ASSERT(suplex_read_lp(solver, fixture("tiny.lp").c_str()) == SUPLEX_OK);
    TEST_ASSERT(suplex_set_algorithm(solver, 1) == SUPLEX_OK);
    TEST_ASSERT(suplex_set_presolve(solver, 0) == SUPLEX_OK);
    TEST_ASSERT(suplex_solve(solver) == static_cast<int>(SolverStatus::OPTIMAL));
    double value = 0.0;
    TEST_ASSERT(suplex_get_solution(solver, &value, 0) == SUPLEX_ERROR_BUFFER_TOO_SMALL);
    suplex_destroy(solver);
}
REGISTER_TEST(test_c_api_validation_and_buffer_size);

void test_solution_writer() {
    Solution solution;
    solution.status = SolverStatus::OPTIMAL;
    solution.objective_value = 7.5;
    solution.primal_values = {2.0};
    solution.row_activities = {2.0};
    solution.dual_values = {1.25};
    const auto path = std::filesystem::temp_directory_path() / "suplex-integration.sol";
    SolutionWriter writer;
    TEST_ASSERT(writer.write(path.string(), solution, {"flow"}, {"balance"}));
    TEST_ASSERT(std::filesystem::file_size(path) > 0);
    std::filesystem::remove(path);
}
REGISTER_TEST(test_solution_writer);

void test_missing_ipm_reports_error() {
    Suplex solver;
    TEST_ASSERT(solver.read_lp(fixture("tiny.lp")));
    solver.set_algorithm(Algorithm::IPM);
    TEST_ASSERT(solver.solve() == SolverStatus::NUMERICAL_ERROR);
    TEST_ASSERT(solver.last_error().find("P2") != std::string::npos);
}
REGISTER_TEST(test_missing_ipm_reports_error);

void test_milp_end_to_end() {
    Suplex solver;
    TEST_ASSERT(solver.read_lp(fixture("tiny_mip.lp")));
    solver.set_algorithm(Algorithm::PRIMAL_SIMPLEX);
    solver.set_log_level(LogLevel::OFF);
    const SolverStatus status = solver.solve();
    TEST_ASSERT(status == SolverStatus::OPTIMAL);
    TEST_ASSERT_NEAR(solver.solution().objective_value, -1.0, 1e-7);
    TEST_ASSERT(solver.solution().primal_values.size() == 2);
}
REGISTER_TEST(test_milp_end_to_end);

void test_milp_fractional_root_branches() {
    Suplex solver;
    TEST_ASSERT(solver.read_lp(fixture("branch_mip.lp")));
    solver.set_algorithm(Algorithm::PRIMAL_SIMPLEX);
    solver.set_log_level(LogLevel::OFF);
    const SolverStatus status = solver.solve();
    TEST_ASSERT(status == SolverStatus::OPTIMAL);
    TEST_ASSERT_NEAR(solver.solution().objective_value, 0.0, 1e-7);
    TEST_ASSERT(solver.solution().primal_values.size() == 1);
    TEST_ASSERT_NEAR(solver.solution().primal_values[0], 0.0, 1e-7);
    TEST_ASSERT(solver.solution().num_nodes >= 1);
}
REGISTER_TEST(test_milp_fractional_root_branches);

void test_gpu_cpu_fallback() {
    const std::vector<Index> rows = {0, 0, 1, 1};
    const std::vector<Index> cols = {0, 1, 0, 1};
    const std::vector<Real> values = {4.0, 1.0, 1.0, 3.0};
    const SparseMatrixCSC matrix = SparseMatrixCSC::from_triplets(2, 2, rows, cols, values);
    const SparseMatrixCSR csr = matrix.to_csr();
    const Real x[] = {1.0, 2.0}; Real y[] = {0.0, 0.0};
    GPUManager::instance().spmv_csr(csr, x, y);
    TEST_ASSERT_NEAR(y[0], 6.0, 1e-12); TEST_ASSERT_NEAR(y[1], 7.0, 1e-12);
    const Real rhs[] = {1.0, 2.0}; Real solution[] = {0.0, 0.0};
    TEST_ASSERT(GPUManager::instance().pcg_solve(matrix, rhs, solution, 1e-12, 100));
    TEST_ASSERT_NEAR(solution[0], 1.0 / 11.0, 1e-10);
    TEST_ASSERT_NEAR(solution[1], 7.0 / 11.0, 1e-10);
}
REGISTER_TEST(test_gpu_cpu_fallback);
} // namespace
