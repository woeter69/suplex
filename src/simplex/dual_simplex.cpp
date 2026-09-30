#include "dual_simplex.h"
#include "primal_simplex.h"
#include <chrono>
#include <cmath>
#include <iostream>
#include <algorithm>

namespace suplex {

DualSimplexSolver::DualSimplexSolver() = default;

void DualSimplexSolver::setup_augmented_system(const Problem& problem) {
    Index m = problem.num_rows();
    Index n = problem.num_cols();
    Index total_vars = n + m;

    TripletMatrix trip(m, total_vars);
    trip.reserve(problem.constraint_matrix().nnz + m);

    const auto& A = problem.constraint_matrix();
    for (Index j = 0; j < n; ++j) {
        Index begin = A.col_begin(j);
        Index end = A.col_end(j);
        for (Index k = begin; k < end; ++k) {
            trip.add_entry(A.row_index[k], j, A.values[k]);
        }
    }

    for (Index i = 0; i < m; ++i) {
        trip.add_entry(i, n + i, -1.0);
    }
    A_aug_ = trip.to_csc();

    c_aug_.assign(total_vars, 0.0);
    const Real sign = (problem.sense() == ObjectiveSense::MAXIMIZE) ? -1.0 : 1.0;
    for (Index j = 0; j < n; ++j) {
        c_aug_[j] = sign * problem.objective()[j];
    }

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

void DualSimplexSolver::compute_basic_primal() {
    Index m = basis_.num_rows();
    Index total_vars = basis_.total_vars();

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

    std::vector<Real> rhs(m, 0.0);
    for (Index j = 0; j < total_vars; ++j) {
        if (!basis_.is_basic(j) && std::abs(x_[j]) > EPS_ZERO) {
            A_aug_.multiply_col(j, -x_[j], rhs.data());
        }
    }

    lu_manager_.ftran(rhs.data(), x_B_.data());
    for (Index i = 0; i < m; ++i) {
        x_[basis_.basic_var_at(i)] = x_B_[i];
    }
}

void DualSimplexSolver::compute_duals_and_reduced_costs() {
    Index m = basis_.num_rows();
    Index total_vars = basis_.total_vars();

    std::vector<Real> c_B(m, 0.0);
    for (Index i = 0; i < m; ++i) {
        c_B[i] = c_aug_[basis_.basic_var_at(i)];
    }

    lu_manager_.btran(c_B.data(), y_.data());

    for (Index j = 0; j < total_vars; ++j) {
        if (basis_.is_basic(j)) {
            d_[j] = 0.0;
        } else {
            d_[j] = c_aug_[j] - A_aug_.dot_col(j, y_.data());
        }
    }
}

void DualSimplexSolver::ensure_dual_feasibility() {
    Index total_vars = basis_.total_vars();
    compute_duals_and_reduced_costs();

    bool flipped = false;
    for (Index j = 0; j < total_vars; ++j) {
        if (basis_.is_basic(j) || basis_.status(j) == BasisStatus::FIXED) continue;

        if (basis_.status(j) == BasisStatus::AT_LOWER && d_[j] < -EPS_OPTIMALITY) {
            if (upper_aug_[j] < INF) {
                basis_.set_status(j, BasisStatus::AT_UPPER);
                x_[j] = upper_aug_[j];
                flipped = true;
            }
        } else if (basis_.status(j) == BasisStatus::AT_UPPER && d_[j] > EPS_OPTIMALITY) {
            if (lower_aug_[j] > -INF) {
                basis_.set_status(j, BasisStatus::AT_LOWER);
                x_[j] = lower_aug_[j];
                flipped = true;
            }
        }
    }

    if (flipped) {
        compute_basic_primal();
    }
}

void DualSimplexSolver::extract_solution(
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

        Real obj = 0.0;
        for (Index j = 0; j < n; ++j) {
            obj += problem.objective()[j] * solution.primal_values[j];
        }
        solution.objective_value = obj;

        solution.row_activities.assign(m, 0.0);
        problem.constraint_matrix().multiply_vec(solution.primal_values.data(), solution.row_activities.data());
    }
}

SolverStatus DualSimplexSolver::solve(const Problem& problem, Solution& solution) {
    auto start_time = std::chrono::steady_clock::now();

    setup_augmented_system(problem);
    Index m = problem.num_rows();
    Index total_vars = problem.num_cols() + m;

    if (basis_.num_rows() != m || basis_.total_vars() != total_vars) {
        basis_.init(m, total_vars);
        basis_.crash(A_aug_, lower_aug_, upper_aug_);
    }

    if (lu_manager_.refactorize(A_aug_, basis_.basic_vars()) == SolverStatus::NUMERICAL_ERROR) {
        basis_.set_logical_basis(problem.num_cols());
        if (lu_manager_.refactorize(A_aug_, basis_.basic_vars()) == SolverStatus::NUMERICAL_ERROR) {
            solution.status = SolverStatus::NUMERICAL_ERROR;
            return SolverStatus::NUMERICAL_ERROR;
        }
    }

    compute_basic_primal();
    ensure_dual_feasibility();

    // If starting from a cold basis that cannot be made dual feasible,
    // solve via PrimalSimplexSolver.
    bool dual_feasible = true;
    for (Index j = 0; j < total_vars; ++j) {
        if (basis_.is_basic(j) || basis_.status(j) == BasisStatus::FIXED) continue;
        if (basis_.status(j) == BasisStatus::AT_LOWER && d_[j] < -EPS_OPTIMALITY) {
            dual_feasible = false;
            break;
        }
        if (basis_.status(j) == BasisStatus::AT_UPPER && d_[j] > EPS_OPTIMALITY) {
            dual_feasible = false;
            break;
        }
    }

    if (!dual_feasible) {
        PrimalSimplexSolver primal_solver;
        primal_solver.set_iteration_limit(max_iter_);
        primal_solver.set_time_limit(time_limit_);
        primal_solver.set_log_level(log_level_);
        SolverStatus stat = primal_solver.solve(problem, solution);
        // Sync basis
        basis_.init(m, total_vars);
        const auto& pbasis = primal_solver.get_basis();
        for (Index i = 0; i < m && i < static_cast<Index>(pbasis.size()); ++i) {
            basis_.pivot(i, pbasis[i], BasisStatus::AT_LOWER);
        }
        return stat;
    }

    int64_t iter = 0;
    SolverStatus final_status = SolverStatus::NOT_STARTED;

    while (iter < max_iter_) {
        auto current_time = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(current_time - start_time).count();
        if (elapsed > time_limit_) {
            final_status = SolverStatus::TIME_LIMIT;
            break;
        }

        // Step 1: Select leaving variable (primal infeasibility)
        Index leaving_row = -1;
        Real max_infeas = EPS_FEASIBILITY;
        Real leaving_dir = 0.0; // +1 if below lower, -1 if above upper

        for (Index i = 0; i < m; ++i) {
            Index var = basis_.basic_var_at(i);
            Real val = x_B_[i];
            Real lb = lower_aug_[var];
            Real ub = upper_aug_[var];

            Real below = lb - val;
            Real above = val - ub;

            if (below > max_infeas) {
                max_infeas = below;
                leaving_row = i;
                leaving_dir = 1.0; // needs to increase
            } else if (above > max_infeas) {
                max_infeas = above;
                leaving_row = i;
                leaving_dir = -1.0; // needs to decrease
            }
        }

        if (leaving_row == -1) {
            // No primal infeasibilities remain -> OPTIMAL
            final_status = SolverStatus::OPTIMAL;
            break;
        }

        // Step 2: Compute pivot row
        // tau = B^{-T} * e_{leaving_row}
        std::vector<Real> e_p(m, 0.0);
        e_p[leaving_row] = 1.0;
        std::vector<Real> tau(m, 0.0);
        lu_manager_.btran(e_p.data(), tau.data());

        // Compute row nonzeros: rho_j = tau^T * A_aug[:, j]
        std::vector<Real> rho(total_vars, 0.0);
        for (Index j = 0; j < total_vars; ++j) {
            if (!basis_.is_basic(j)) {
                rho[j] = A_aug_.dot_col(j, tau.data());
            }
        }

        // Step 3: Dual ratio test (Harris)
        Index best_q = -1;
        Real harris_theta = INF;
        const Real HARRIS_TOL = 1e-7;

        for (Index j = 0; j < total_vars; ++j) {
            if (basis_.is_basic(j) || basis_.status(j) == BasisStatus::FIXED) continue;

            BasisStatus st = basis_.status(j);
            Real dj = d_[j];
            Real delta_j = -leaving_dir * rho[j];

            if (st == BasisStatus::AT_LOWER && delta_j > EPS_PIVOT) {
                Real ratio = (std::max<Real>(0.0, dj) + HARRIS_TOL) / delta_j;
                harris_theta = std::min(harris_theta, ratio);
            } else if (st == BasisStatus::AT_UPPER && delta_j < -EPS_PIVOT) {
                Real ratio = (std::max<Real>(0.0, -dj) + HARRIS_TOL) / (-delta_j);
                harris_theta = std::min(harris_theta, ratio);
            }
        }

        if (harris_theta >= INF / 2.0) {
            // Dual unbounded => Primal infeasible
            final_status = SolverStatus::INFEASIBLE;
            break;
        }

        // Pass 2: Select candidate with largest |delta_j| within harris_theta
        Real best_pivot_abs = 0.0;
        for (Index j = 0; j < total_vars; ++j) {
            if (basis_.is_basic(j) || basis_.status(j) == BasisStatus::FIXED) continue;

            BasisStatus st = basis_.status(j);
            Real dj = d_[j];
            Real delta_j = -leaving_dir * rho[j];

            if (st == BasisStatus::AT_LOWER && delta_j > EPS_PIVOT) {
                Real ratio = std::max<Real>(0.0, dj) / delta_j;
                if (ratio <= harris_theta && delta_j > best_pivot_abs) {
                    best_pivot_abs = delta_j;
                    best_q = j;
                }
            } else if (st == BasisStatus::AT_UPPER && delta_j < -EPS_PIVOT) {
                Real ratio = std::max<Real>(0.0, -dj) / (-delta_j);
                if (ratio <= harris_theta && (-delta_j) > best_pivot_abs) {
                    best_pivot_abs = -delta_j;
                    best_q = j;
                }
            }
        }

        if (best_q == -1) {
            final_status = SolverStatus::INFEASIBLE;
            break;
        }

        // Step 4: Compute pivot column
        std::vector<Real> a_q(m, 0.0);
        Index begin = A_aug_.col_begin(best_q);
        Index end = A_aug_.col_end(best_q);
        for (Index k = begin; k < end; ++k) {
            a_q[A_aug_.row_index[k]] = A_aug_.values[k];
        }

        std::vector<Real> alpha(m, 0.0);
        lu_manager_.ftran(a_q.data(), alpha.data());

        Real pivot_elem = alpha[leaving_row];
        if (std::abs(pivot_elem) < EPS_PIVOT) {
            final_status = SolverStatus::NUMERICAL_ERROR;
            break;
        }

        // Step 5: Update
        Index leaving_var = basis_.basic_var_at(leaving_row);
        BasisStatus leaving_status = (leaving_dir > 0) ? BasisStatus::AT_LOWER : BasisStatus::AT_UPPER;
        x_[leaving_var] = (leaving_status == BasisStatus::AT_LOWER) ? lower_aug_[leaving_var] : upper_aug_[leaving_var];

        // Perform basis pivot
        basis_.pivot(leaving_row, best_q, leaving_status);

        // Update LU
        SolverStatus lu_stat = lu_manager_.update(leaving_row, alpha);
        if (lu_stat == SolverStatus::NUMERICAL_ERROR || lu_manager_.should_refactorize()) {
            if (lu_manager_.refactorize(A_aug_, basis_.basic_vars()) == SolverStatus::NUMERICAL_ERROR) {
                final_status = SolverStatus::NUMERICAL_ERROR;
                break;
            }
        }

        compute_basic_primal();
        compute_duals_and_reduced_costs();

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

SolverStatus DualSimplexSolver::solve_from_basis(
    const Problem& problem,
    const std::vector<Index>& basis_indices,
    Solution& solution
) {
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
