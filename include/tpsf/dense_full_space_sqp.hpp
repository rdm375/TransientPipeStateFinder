#pragma once

#include "tpsf/dense_constrained_least_squares.hpp"
#include "tpsf/dense_full_space_subproblem.hpp"

#include <Eigen/Dense>

#include <cstddef>
#include <functional>
#include <vector>

namespace tpsf {

enum class DenseSQPTermination {
    converged,
    maximum_iterations,
    incompatible_linearization,
    non_descent_direction,
    line_search_failure,
    numerical_failure
};

struct DenseSQPOptions {
    std::size_t maximum_iterations = 40;
    std::size_t maximum_backtracks = 24;

    double feasibility_tolerance = 1.0e-8;
    double stationarity_tolerance = 1.0e-8;

    double initial_penalty = 1.0;
    double penalty_margin = 0.1;

    double backtracking_factor = 0.5;
    double armijo_fraction = 1.0e-4;

    DenseConstrainedLeastSquaresOptions linear_solver;
};

struct DenseSQPIteration {
    std::size_t iteration = 0;

    double objective = 0.0;
    double constraint_l1 = 0.0;
    double constraint_l2 = 0.0;
    double projected_gradient = 0.0;

    double step_norm = 0.0;
    double step_length = 0.0;
    double penalty = 0.0;

    std::size_t constraint_rank = 0;
    std::size_t reduced_measurement_rank = 0;
    std::size_t unresolved_directions = 0;
};

struct DenseSQPResult {
    Eigen::VectorXd coordinates;

    DenseSQPTermination termination =
        DenseSQPTermination::numerical_failure;

    std::vector<DenseSQPIteration> history;

    double final_objective = 0.0;
    double final_constraint_l2 = 0.0;
    double final_projected_gradient = 0.0;

    bool converged() const noexcept {
        return termination == DenseSQPTermination::converged;
    }
};

using DenseSQPAssembly =
    std::function<DenseFullSpaceSubproblem(
        const Eigen::VectorXd&)>;

DenseSQPResult solve_dense_full_space_sqp(
    const DenseSQPAssembly& assemble,
    const Eigen::VectorXd& initial,
    const DenseSQPOptions& options = {});

DenseSQPResult solve_dense_full_space_sqp(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& initial,
    const DenseSQPOptions& options = {});

} // namespace tpsf
