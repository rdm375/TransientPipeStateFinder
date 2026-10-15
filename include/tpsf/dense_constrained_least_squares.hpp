#pragma once

#include <Eigen/Dense>

#include <cstddef>

namespace tpsf {

struct DenseConstrainedLeastSquaresOptions {
    // Relative singular-value cutoff:
    // sigma_i > relative_rank_tolerance * sigma_max.
    //
    // Zero singular values are always excluded.
    double relative_rank_tolerance = 1.0e-12;

    // Relative compatibility tolerance for equality constraints.
    double feasibility_tolerance = 1.0e-10;
};

struct DenseConstrainedLeastSquaresResult {
    Eigen::VectorXd step;

    std::size_t constraint_rank = 0;
    std::size_t reduced_measurement_rank = 0;
    std::size_t unresolved_directions = 0;

    bool linearized_constraints_compatible = false;

    double linearized_constraint_residual = 0.0;
    double linearized_measurement_residual = 0.0;

    double constraint_singular_value_max = 0.0;
    double constraint_singular_value_min_retained = 0.0;

    double reduced_measurement_singular_value_max = 0.0;
    double reduced_measurement_singular_value_min_retained = 0.0;
};

DenseConstrainedLeastSquaresResult
solve_dense_constrained_least_squares(
    const Eigen::MatrixXd& constraint_jacobian,
    const Eigen::VectorXd& constraint_residual,
    const Eigen::MatrixXd& measurement_jacobian,
    const Eigen::VectorXd& measurement_residual,
    const DenseConstrainedLeastSquaresOptions& options = {});

} // namespace tpsf
