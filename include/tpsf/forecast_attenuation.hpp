#pragma once

#include "tpsf/forecast_continuation.hpp"
#include "tpsf/full_space_problem.hpp"
#include "tpsf/statefinding_quality.hpp"

#include <Eigen/Dense>
#include <vector>

namespace tpsf {

struct ForecastSensitivityPoint {
    double time = 0.0;
    double linepack = 0.0;
    double outlet_pressure = 0.0;

    // Derivatives per unit scaled reconstruction coordinate.
    double linepack_feasible = 0.0;
    double linepack_unresolved = 0.0;
    double outlet_pressure_feasible = 0.0;
    double outlet_pressure_unresolved = 0.0;
};

struct ForecastAttenuationResult {
    std::vector<ForecastSensitivityPoint> points;
    double maximum_mass_defect = 0.0;
    double maximum_tangent_residual = 0.0;
};

// One nominal nonlinear forecast; all directional derivatives
// use exact discrete JVPs and implicit linear solves.
//
// Future boundary inputs are held fixed.
ForecastAttenuationResult assess_forecast_attenuation(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& coordinates,
    const StatefindingQuality& quality,
    const Configuration& forecast_configuration,
    const BoundarySchedule& future_boundary);

} // namespace tpsf
