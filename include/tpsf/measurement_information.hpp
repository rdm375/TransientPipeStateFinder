#pragma once

#include <Eigen/Dense>

#include <cstddef>

namespace tpsf {

struct MeasurementInformationDecomposition {
    Eigen::MatrixXd resolved_basis;
    Eigen::MatrixXd unresolved_basis;
    Eigen::VectorXd retained_singular_values;

    std::size_t feasible_dimension = 0;
    std::size_t information_rank = 0;
    std::size_t unresolved_dimension = 0;
};

// Decompose standardized measurement information restricted
// to the physically feasible tangent subspace.
//
// B: measurement Jacobian in scaled optimization coordinates.
// Z: orthonormal basis of the constraint-feasible subspace.
//
// The returned bases are in scaled optimization coordinates.
MeasurementInformationDecomposition
decompose_measurement_information(
    const Eigen::MatrixXd& B,
    const Eigen::MatrixXd& Z,
    double relative_rank_tolerance = 1.0e-12);

} // namespace tpsf
