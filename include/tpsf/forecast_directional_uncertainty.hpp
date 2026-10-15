#pragma once

#include "tpsf/forecast_attenuation.hpp"
#include "tpsf/measurement_information.hpp"

#include <Eigen/Dense>
#include <vector>

namespace tpsf {

struct ForecastDirectionalUncertainty {
    std::vector<double> time;
    std::vector<double> nominal_outlet_pressure;

    // Rows correspond to forecast times.
    // Columns correspond to singular-vector directions.
    Eigen::MatrixXd resolved_sensitivity;
    Eigen::MatrixXd unresolved_sensitivity;

    Eigen::VectorXd retained_singular_values;

    // Radius for unit standardized measurement perturbation.
    std::vector<double> measurement_unit_radius;

    // Radius for unit unresolved-coordinate perturbation.
    std::vector<double> unresolved_unit_radius;

    // True when the linearized forecast has nonzero
    // sensitivity to an unresolved direction.
    //
    // This is an algebraic diagnostic, not a numerical
    // significance test or nonlinear certificate.
    std::vector<bool> unresolved_influence;

    double maximum_tangent_residual = 0.0;
    double maximum_tangent_mass_defect = 0.0;
    double maximum_nominal_mass_defect = 0.0;
};

ForecastDirectionalUncertainty
assess_forecast_directional_uncertainty(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& coordinates,
    const MeasurementInformationDecomposition& information,
    const Configuration& configuration,
    const BoundarySchedule& boundary);

} // namespace tpsf
