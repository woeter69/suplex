#pragma once

#include "lp_solver.h"
#include "basis.h"
#include "lu_update.h"
#include "pricing.h"

namespace suplex {

/**
 * @brief Revised Primal Simplex Solver for general bounded linear programs.
 */
class PrimalSimplexSolver : public LPSolver {
public:
    PrimalSimplexSolver();
    ~PrimalSimplexSolver() override = default;

    SolverStatus solve(const Problem& problem, Solution& solution) override;

    SolverStatus solve_from_basis(
        const Problem& problem,
        const std::vector<Index>& basis_indices,
        Solution& solution
    ) override;

    void set_iteration_limit(Index limit) override { max_iter_ = limit; }
    void set_time_limit(Real seconds) override { time_limit_ = seconds; }
    void set_log_level(LogLevel level) override { log_level_ = level; }
    void set_pricing_strategy(PricingStrategy strat) { pricing_strategy_ = strat; }

    const std::vector<Index>& get_basis() const override { return basis_.basic_vars(); }

private:
    Index max_iter_ = MAX_ITER;
    Real time_limit_ = TIME_LIMIT;
    LogLevel log_level_ = LogLevel::INFO;
    PricingStrategy pricing_strategy_ = PricingStrategy::DEVEX;

    Basis basis_;
    LUManager lu_manager_;
    Pricing pricing_;

    // Internal data structures for augmented system [A, -I]
    SparseMatrixCSC A_aug_;
    std::vector<Real> c_aug_;
    std::vector<Real> lower_aug_;
    std::vector<Real> upper_aug_;
    std::vector<Real> x_;
    std::vector<Real> x_B_;
    std::vector<Real> y_;
    std::vector<Real> d_;

    void setup_augmented_system(const Problem& problem);
    void compute_basic_solution();
    void compute_duals_and_reduced_costs(const std::vector<Real>& cost_vec);
    bool check_primal_feasibility(Real& max_infeas) const;
    void extract_solution(const Problem& problem, Solution& solution, SolverStatus status, int64_t iters, double elapsed);
};

} // namespace suplex
