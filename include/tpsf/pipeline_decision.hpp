#pragma once

#include "tpsf/decision_reliability.hpp"
#include "tpsf/forecast_attenuation.hpp"
#include "tpsf/forecast_directional_uncertainty.hpp"

namespace tpsf {

// Conservative local-linearized uncertainty model.
//
// measurement_radius is a bound on the Euclidean norm of
// standardized measurement perturbations. It is NOT
// automatically a statistical confidence level.
//
// If bound_unresolved is false, unresolved directions are
// unbounded. Otherwise unresolved_radius bounds their
// Euclidean displacement in scaled coordinates.
struct PipelineDecisionOptions {
    double minimum_outlet_pressure_pa = 0.0;
    double measurement_radius = 1.0;

    bool bound_unresolved = false;
    double unresolved_radius = 0.0;
};

struct PipelineDecisionResult {
    ForecastAttenuationResult attenuation;
    ForecastDirectionalUncertainty directional;
    ForecastDecisionResult decision;

    // Conservative displacement radius in the feasible
    // coordinate subspace induced by measurement uncertainty.
    double measurement_coordinate_radius = 0.0;

    bool unresolved_bounded = false;
};

// Assess minimum outlet pressure at sampled forecast times.
// Future boundary inputs are fixed.
//
// The uncertainty bounds are conservative for the
// linearized model; they are not nonlinear certificates.
PipelineDecisionResult assess_pipeline_decision(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& coordinates,
    const StatefindingQuality& quality,
    const Configuration& forecast_configuration,
    const BoundarySchedule& future_boundary,
    const PipelineDecisionOptions& options);

} // namespace tpsf
