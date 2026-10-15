#pragma once

#include "tpsf/engineering_sensitivity.hpp"

#include <Eigen/Dense>
#include <cstddef>

namespace tpsf {

enum class UncertaintyStatus {
    finite,
    unbounded,
    nondifferentiable
};

struct MeasurementUncertaintyOptions {
    // Relative SVD rank tolerance.
    double relative_rank_tolerance = 1.0e-12;

    // Relative threshold for detecting engineering sensitivity
    // along unresolved directions.
    double unresolved_relative_tolerance = 1.0e-10;

    // Absolute sensitivity threshold in engineering units
    // per unit scaled-coordinate displacement.
    double unresolved_absolute_tolerance = 1.0e-12;
};

struct EngineeringUncertainty {
    UncertaintyStatus status = UncertaintyStatus::unbounded;

    // Local, linearized, measurement-noise-only uncertainty.
    // Valid as a complete data-only result only when finite.
    double variance = 0.0;
    double standard_deviation = 0.0;

    // Component of engineering sensitivity not resolved
    // by the measurement data.
    double unresolved_sensitivity_norm = 0.0;

    std::size_t feasible_dimension = 0;
    std::size_t information_rank = 0;
    std::size_t unresolved_directions = 0;

    // False when an engineering quantity depends on an
    // unresolved direction.
    bool finite_data_only_bound = false;
};

// General analytic interface.
//
// S: standardized measurement Jacobian restricted to
//    feasible coordinates, with unit-variance residual noise.
//
// g: engineering-quantity gradient in the same feasible
//    coordinates.
//
// The function detects rank deficiency and never silently
// interprets a pseudoinverse as a complete uncertainty bound.
EngineeringUncertainty propagate_measurement_uncertainty(
    const Eigen::MatrixXd& information_jacobian,
    const Eigen::VectorXd& feasible_quantity_gradient,
    const MeasurementUncertaintyOptions& options = {});

// Adapter for E5.1 and E5.2.
//
// The problem is required because the standardized measurement
// Jacobian must be assembled at the specified coordinates.
EngineeringUncertainty assess_engineering_uncertainty(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& coordinates,
    const StatefindingQuality& quality,
    const QuantitySensitivity& sensitivity,
    const MeasurementUncertaintyOptions& options = {});

} // namespace tpsf
