#pragma once

#include "tpsf/full_space_problem.hpp"

#include <Eigen/Dense>

namespace tpsf {

struct DenseFullSpaceSubproblem {
    Eigen::MatrixXd constraint_jacobian;
    Eigen::VectorXd constraint_residual;

    Eigen::MatrixXd measurement_jacobian;
    Eigen::VectorXd measurement_residual;
};

// Assemble the exact linearized, scaled constrained
// least-squares subproblem at optimization coordinates u.
//
// No finite differences are used.
DenseFullSpaceSubproblem assemble_dense_full_space_subproblem(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& u);

} // namespace tpsf
