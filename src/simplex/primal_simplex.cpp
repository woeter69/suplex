#include "primal_simplex.h"
#include <chrono>
#include <cmath>
#include <iostream>
#include <algorithm>

namespace suplex {

PrimalSimplexSolver::PrimalSimplexSolver() = default;

void PrimalSimplexSolver::setup_augmented_system(const Problem& problem) {
    Index m = problem.num_rows();
    Index n = problem.num_cols();
    Index total_vars = n + m;

    TripletMatrix trip(m, total_vars);
    trip.reserve(problem.constraint_matrix().nnz + m);

    // Copy structural columns
    const auto& A = problem.constraint_matrix();
    for (Index j = 0; j < n; ++j) {
        Index begin = A.col_begin(j);
        Index end = A.col_end(j);
        for (Index k = begin; k < end; ++k) {
            trip.add_entry(A.row_index[k], j, A.values[k]);
        }
    }

    // Add logical slack columns: -I_m
    for (Index i = 0; i < m; ++i) {
        trip.add_entry(i, n + i, -1.0);
    }
    A_aug_ = trip.to_csc();

    // Objective
    c_aug_.assign(total_vars, 0.0);
    const Real sign = (problem.sense() == ObjectiveSense::MAXIMIZE) ? -1.0 : 1.0;
    for (Index j = 0; j < n; ++j) {
        c_aug_[j] = sign * problem.objective()[j];
    }

    // Bounds
    lower_aug_.resize(total_vars);
    upper_aug_.resize(total_vars);
    for (Index j = 0; j < n; ++j) {
        lower_aug_[j] = problem.col_lower()[j];
        upper_aug_[j] = problem.col_upper()[j];
    }
    for (Index i = 0; i < m; ++i) {
        lower_aug_[n + i] = problem.row_lower()[i];
        upper_aug_[n + i] = problem.row_upper()[i];
    }

    x_.assign(total_vars, 0.0);
    x_B_.assign(m, 0.0);
    y_.assign(m, 0.0);
    d_.assign(total_vars, 0.0);
}

void PrimalSimplexSolver::compute_basic_solution() {
    Index m = basis_.num_rows();
    Index total_vars = basis_.total_vars();

    // Set non-basic variables to their bound
    for (Index j = 0; j < total_vars; ++j) {
        BasisStatus st = basis_.status(j);
        if (st == BasisStatus::BASIC) continue;

        if (st == BasisStatus::AT_LOWER) {
            x_[j] = (lower_aug_[j] > -INF) ? lower_aug_[j] : 0.0;
        } else if (st == BasisStatus::AT_UPPER) {
            x_[j] = (upper_aug_[j] < INF) ? upper_aug_[j] : 0.0;
        } else if (st == BasisStatus::FIXED) {
            x_[j] = lower_aug_[j];
        } else {
            x_[j] = 0.0;
        }
    }

    // rhs = - sum_{j in N} A_aug[:, j] * x_[j]
    std::vector<Real> rhs(m, 0.0);
    for (Index j = 0; j < total_vars; ++j) {
        if (!basis_.is_basic(j) && std::abs(x_[j]) > EPS_ZERO) {
            A_aug_.multiply_col(j, -x_[j], rhs.data());
        }
    }

    // Solve B * x_B = rhs
    lu_manager_.ftran(rhs.data(), x_B_.data());
    for (Index i = 0; i < m; ++i) {
        x_[basis_.basic_var_at(i)] = x_B_[i];
    }
}

void PrimalSimplexSolver::compute_duals_and_reduced_costs(const std::vector<Real>& cost_vec) {
    Index m = basis_.num_rows();
    Index total_vars = basis_.total_vars();

    std::vector<Real> c_B(m, 0.0);
    for (Index i = 0; i < m; ++i) {
        c_B[i] = cost_vec[basis_.basic_var_at(i)];
    }

    // Solve B^T * y = c_B
    lu_manager_.btran(c_B.data(), y_.data());

    // Compute reduced costs: d_j = c_j - y^T A_aug[:, j]
    for (Index j = 0; j < total_vars; ++j) {
        if (basis_.is_basic(j)) {
            d_[j] = 0.0;
        } else {
            d_[j] = cost_vec[j] - A_aug_.dot_col(j, y_.data());
        }
    }
}

bool PrimalSimplexSolver::check_primal_feasibility(Real& max_infeas) const {
    Index m = basis_.num_rows();
    max_infeas = 0.0;
    bool feasible = true;

    for (Index i = 0; i < m; ++i) {
        Index var = basis_.basic_var_at(i);
        Real val = x_B_[i];
        Real lb = lower_aug_[var];
        Real ub = upper_aug_[var];

        if (val < lb - EPS_FEASIBILITY) {
            feasible = false;
            max_infeas = std::max(max_infeas, lb - val);
        } else if (val > ub + EPS_FEASIBILITY) {
            feasible = false;
            max_infeas = std::max(max_infeas, val - ub);
        }
    }
    return feasible;
}

void PrimalSimplexSolver::extract_solution(
    const Problem& problem,
    Solution& solution,
    SolverStatus status,
    int64_t iters,
    double elapsed
) {
    Index n = problem.num_cols();
    Index m = problem.num_rows();

    solution.status = status;
    solution.num_iterations = iters;
    solution.solve_time_seconds = elapsed;

    if (status == SolverStatus::OPTIMAL || status == SolverStatus::ITERATION_LIMIT || status == SolverStatus::TIME_LIMIT) {
        solution.primal_values.assign(x_.begin(), x_.begin() + n);
        solution.dual_values = y_;
        solution.reduced_costs.assign(d_.begin(), d_.begin() + n);

        const Real sign = (problem.sense() == ObjectiveSense::MAXIMIZE) ? -1.0 : 1.0;
        for (Index i = 0; i < m; ++i) {
            solution.dual_values[i] *= sign;
        }
        for (Index j = 0; j < n; ++j) {
            solution.reduced_costs[j] *= sign;
        }

        // Calculate original objective value
        Real obj = 0.0;
        for (Index j = 0; j < n; ++j) {
            obj += problem.objective()[j] * solution.primal_values[j];
        }
        solution.objective_value = obj;

        // Calculate row activities Ax
        solution.row_activities.assign(m, 0.0);
        problem.constraint_matrix().multiply_vec(solution.primal_values.data(), solution.row_activities.data());
    }
}

SolverStatus PrimalSimplexSolver::solve(const Problem& problem, Solution& solution) {
    auto start_time = std::chrono::steady_clock::now();

    setup_augmented_system(problem);
    Index m = problem.num_rows();
    Index total_vars = problem.num_cols() + m;

    basis_.init(m, total_vars);
    basis_.crash(A_aug_, lower_aug_, upper_aug_);

    pricing_.init(total_vars, m, pricing_strategy_);

    if (lu_manager_.refactorize(A_aug_, basis_.basic_vars()) == SolverStatus::NUMERICAL_ERROR) {
        // Fall back to pure logical basis
        basis_.set_logical_basis(problem.num_cols());
        if (lu_manager_.refactorize(A_aug_, basis_.basic_vars()) == SolverStatus::NUMERICAL_ERROR) {
            solution.status = SolverStatus::NUMERICAL_ERROR;
            return SolverStatus::NUMERICAL_ERROR;
        }
    }

    compute_basic_solution();

    Real max_infeas = 0.0;
    bool phase1 = !check_primal_feasibility(max_infeas);

    int64_t iter = 0;
    int degenerate_count = 0;
    SolverStatus final_status = SolverStatus::NOT_STARTED;

    while (iter < max_iter_) {
        auto current_time = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(current_time - start_time).count();
        if (elapsed > time_limit_) {
            final_status = SolverStatus::TIME_LIMIT;
            break;
        }

        // Step 1: Objective vector for current phase
        std::vector<Real> current_cost(total_vars, 0.0);
        if (phase1) {
            for (Index i = 0; i < m; ++i) {
                Index var = basis_.basic_var_at(i);
                if (x_B_[i] < lower_aug_[var] - EPS_FEASIBILITY) {
                    current_cost[var] = -1.0;
                } else if (x_B_[i] > upper_aug_[var] + EPS_FEASIBILITY) {
                    current_cost[var] = 1.0;
                }
            }
        } else {
            current_cost = c_aug_;
        }

        // Step 2: Compute duals and reduced costs
        compute_duals_and_reduced_costs(current_cost);

        // Step 3: Pricing
        Index q = pricing_.select_entering(d_, basis_.var_status());

        if (q == -1) {
            // No eligible entering variable found
            if (phase1) {
                Real remaining_infeas = 0.0;
                if (!check_primal_feasibility(remaining_infeas) && remaining_infeas > EPS_FEASIBILITY) {
                    final_status = SolverStatus::INFEASIBLE;
                    break;
                }
                // Transition to Phase 2
                phase1 = false;
                pricing_.reset_weights();
                continue;
            } else {
                final_status = SolverStatus::OPTIMAL;
                break;
            }
        }

        // Step 4: Compute pivot column alpha = B^{-1} A_aug[:, q]
        std::vector<Real> a_q(m, 0.0);
        Index begin = A_aug_.col_begin(q);
        Index end = A_aug_.col_end(q);
        for (Index k = begin; k < end; ++k) {
            a_q[A_aug_.row_index[k]] = A_aug_.values[k];
        }

        std::vector<Real> alpha(m, 0.0);
        lu_manager_.ftran(a_q.data(), alpha.data());

        // Step 5: Ratio test
        Real direction = 1.0;
        BasisStatus st = basis_.status(q);
        if (st == BasisStatus::AT_UPPER) {
            direction = -1.0;
            for (Index i = 0; i < m; ++i) alpha[i] = -alpha[i];
        } else if (st == BasisStatus::FREE_ZERO && d_[q] > 0.0) {
            direction = -1.0;
            for (Index i = 0; i < m; ++i) alpha[i] = -alpha[i];
        }

        Real theta_max = INF;
        Index leaving_row = -1;
        BasisStatus leaving_status = BasisStatus::AT_LOWER;

        for (Index i = 0; i < m; ++i) {
            Index var = basis_.basic_var_at(i);
            Real ai = alpha[i];

            if (ai > EPS_PIVOT) {
                Real slack = x_B_[i] - lower_aug_[var];
                Real ratio = slack / ai;
                if (ratio < theta_max) {
                    theta_max = ratio;
                    leaving_row = i;
                    leaving_status = BasisStatus::AT_LOWER;
                }
            } else if (ai < -EPS_PIVOT) {
                Real slack = upper_aug_[var] - x_B_[i];
                Real ratio = slack / (-ai);
                if (ratio < theta_max) {
                    theta_max = ratio;
                    leaving_row = i;
                    leaving_status = BasisStatus::AT_UPPER;
                }
            }
        }

        // Bound flipping: check entering variable's opposite bound
        Real entering_range = upper_aug_[q] - lower_aug_[q];
        if (entering_range < theta_max) {
            // Flip non-basic variable to opposite bound
            if (st == BasisStatus::AT_LOWER) {
                x_[q] = upper_aug_[q];
                basis_.set_status(q, BasisStatus::AT_UPPER);
            } else {
                x_[q] = lower_aug_[q];
                basis_.set_status(q, BasisStatus::AT_LOWER);
            }
            for (Index i = 0; i < m; ++i) {
                x_B_[i] -= entering_range * alpha[i];
                x_[basis_.basic_var_at(i)] = x_B_[i];
            }
            iter++;
            continue;
        }

        if (leaving_row == -1) {
            if (phase1) {
                final_status = SolverStatus::INFEASIBLE;
            } else {
                final_status = SolverStatus::UNBOUNDED;
            }
            break;
        }

        // Step 6: Execute pivot
        theta_max = std::max<Real>(0.0, theta_max);
        if (theta_max < EPS_ZERO) {
            degenerate_count++;
            if (degenerate_count > 100) {
                pricing_.set_strategy(PricingStrategy::BLAND);
            }
        } else {
            degenerate_count = 0;
        }

        Real new_val = (direction > 0) ? (lower_aug_[q] + theta_max) : (upper_aug_[q] - theta_max);
        x_[q] = new_val;

        for (Index i = 0; i < m; ++i) {
            x_B_[i] -= theta_max * alpha[i];
            x_[basis_.basic_var_at(i)] = x_B_[i];
        }

        Index leaving_var = basis_.basic_var_at(leaving_row);
        x_[leaving_var] = (leaving_status == BasisStatus::AT_LOWER) ? lower_aug_[leaving_var] : upper_aug_[leaving_var];
        x_B_[leaving_row] = new_val;

        basis_.pivot(leaving_row, q, leaving_status);

        // Update pricing weights
        pricing_.update(q, leaving_row, alpha);

        // Step 7: LU update
        SolverStatus lu_stat = lu_manager_.update(leaving_row, alpha);
        if (lu_stat == SolverStatus::NUMERICAL_ERROR || lu_manager_.should_refactorize()) {
            if (lu_manager_.refactorize(A_aug_, basis_.basic_vars()) == SolverStatus::NUMERICAL_ERROR) {
                final_status = SolverStatus::NUMERICAL_ERROR;
                break;
            }
            compute_basic_solution();
        }

        if (phase1) {
            Real infeas = 0.0;
            if (check_primal_feasibility(infeas)) {
                phase1 = false;
                pricing_.reset_weights();
            }
        }

        iter++;
    }

    if (final_status == SolverStatus::NOT_STARTED) {
        final_status = SolverStatus::ITERATION_LIMIT;
    }

    auto end_time = std::chrono::steady_clock::now();
    double total_elapsed = std::chrono::duration<double>(end_time - start_time).count();
    extract_solution(problem, solution, final_status, iter, total_elapsed);
    return final_status;
}

SolverStatus PrimalSimplexSolver::solve_from_basis(
    const Problem& problem,
    const std::vector<Index>& basis_indices,
    Solution& solution
) {
    // If warm basis provided, use it
    setup_augmented_system(problem);
    Index m = problem.num_rows();
    Index total_vars = problem.num_cols() + m;

    basis_.init(m, total_vars);
    for (Index i = 0; i < m; ++i) {
        Index var = basis_indices[i];
        basis_.pivot(i, var, BasisStatus::AT_LOWER);
    }
    return solve(problem, solution);
}

} // namespace suplex
