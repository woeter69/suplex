#pragma once

#include "src/core/problem.h"
#include "src/core/solution.h"

#include <string>
#include <vector>

namespace suplex {

/** Algorithm selection for the top-level solver. */
enum class Algorithm : uint8_t { AUTO = 0, PRIMAL_SIMPLEX = 1, DUAL_SIMPLEX = 2, IPM = 3 };

/** High-level C++ facade that owns a problem and dispatches to the solver engines. */
class Suplex {
public:
    bool read_mps(const std::string& filename);
    bool read_lp(const std::string& filename);
    void set_problem(Problem problem);
    SolverStatus solve();

    const Problem& problem() const { return problem_; }
    const Solution& solution() const { return solution_; }
    const std::vector<std::string>& column_names() const { return column_names_; }
    const std::vector<std::string>& row_names() const { return row_names_; }
    const std::string& last_error() const { return error_; }

    void set_algorithm(Algorithm value) { algorithm_ = value; }
    void set_presolve(bool value) { presolve_ = value; }
    void set_gpu(bool value) { gpu_ = value; }
    void set_log_level(LogLevel value) { log_level_ = value; }
    void set_time_limit(Real value) { time_limit_ = value; }
    void set_threads(int value) { threads_ = value; }
    void set_mip_gap(Real value) { mip_gap_ = value; }

private:
    Problem problem_;
    Solution solution_;
    std::vector<std::string> column_names_;
    std::vector<std::string> row_names_;
    std::string error_;
    Algorithm algorithm_ = Algorithm::AUTO;
    bool presolve_ = true;
    bool gpu_ = false;
    LogLevel log_level_ = LogLevel::INFO;
    Real time_limit_ = TIME_LIMIT;
    Real mip_gap_ = 1e-4;
    int threads_ = 0;
};

} // namespace suplex
