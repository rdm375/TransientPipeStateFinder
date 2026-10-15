#include "tpsf/forecast_directional_uncertainty.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace tpsf {
namespace {

std::vector<double> as_std(const Eigen::VectorXd& x)
{
    return {x.data(), x.data() + x.size()};
}

Eigen::VectorXd as_eigen(const std::vector<double>& x)
{
    return Eigen::Map<const Eigen::VectorXd>(
        x.data(), static_cast<Eigen::Index>(x.size()));
}

void propagate_basis(
    const FullSpaceProblem& problem,
    const FullSpaceTrajectory& nominal,
    const Eigen::MatrixXd& basis,
    const pipe_sim::TransientParameters& parameters,
    Eigen::MatrixXd& sensitivities,
    double& maximum_residual,
    double& maximum_mass_defect)
{
    const Eigen::Index width = 2 * nominal.n + 2;

    for (Eigen::Index j = 0; j < basis.cols(); ++j) {
        const auto direction =
            problem.trajectory_direction(
                as_std(basis.col(j)));

        const Eigen::VectorXd initial =
            as_eigen(direction.states.back());

        if (initial.size() != width)
            throw std::runtime_error(
                "Incorrect terminal tangent dimension");

        // Future inlet pressure is prescribed and fixed.
        const double tolerance =
            1.0e-10 * std::max(1.0, initial.norm());

        if (std::abs(initial[0]) > tolerance)
            throw std::invalid_argument(
                "Direction changes prescribed forecast inlet pressure");

        const auto history = propagate_forecast_tangent(
            nominal, initial, parameters);

        if (history.states.size() !=
            static_cast<std::size_t>(sensitivities.rows()))
            throw std::runtime_error(
                "Incomplete directional forecast");

        maximum_residual = std::max(
            maximum_residual,
            history.maximum_equation_residual);

        maximum_mass_defect = std::max(
            maximum_mass_defect,
            history.maximum_mass_balance_defect);

        for (Eigen::Index k = 0;
             k < sensitivities.rows(); ++k) {
            sensitivities(k, j) =
                history.states[static_cast<std::size_t>(k)]
                    [width - 1];
        }
    }
}

} // namespace

ForecastDirectionalUncertainty
assess_forecast_directional_uncertainty(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& coordinates,
    const MeasurementInformationDecomposition& information,
    const Configuration& configuration,
    const BoundarySchedule& boundary)
{
    const Eigen::Index n =
        static_cast<Eigen::Index>(
            problem.free_variable_count());

    if (coordinates.size() != n ||
        !coordinates.allFinite() ||
        information.resolved_basis.rows() != n ||
        information.unresolved_basis.rows() != n ||
        !information.resolved_basis.allFinite() ||
        !information.unresolved_basis.allFinite())
        throw std::invalid_argument(
            "Invalid directional uncertainty dimensions");

    if (configuration.steps < 1)
        throw std::invalid_argument(
            "Forecast requires at least one step");

    const std::size_t count =
        static_cast<std::size_t>(configuration.steps + 1);

    if (boundary.inlet_pressure.size() != count ||
        boundary.outlet_flow.size() != count)
        throw std::invalid_argument(
            "Invalid future boundary schedule");

    const auto reconstruction =
        problem.trajectory(as_std(coordinates));

    const auto recorded = forecast_from_full_space(
        configuration, reconstruction, boundary);

    if (recorded.summary.info != 0 ||
        recorded.summary.completed_steps != configuration.steps)
        throw std::runtime_error(
            "Nominal directional forecast failed");

    FullSpaceTrajectory nominal;
    nominal.n = configuration.cells;

    for (const auto& snapshot :
         recorded.trajectory.snapshots) {
        nominal.states.push_back(snapshot.state);
        nominal.inputs.push_back(snapshot.inputs);
        nominal.reconstructed.push_back(
            snapshot.physical_channels);
    }

    if (nominal.states.size() != count ||
        recorded.summary.outlet_pressure.size() != count)
        throw std::runtime_error(
            "Incomplete nominal directional forecast");

    ForecastDirectionalUncertainty result;

    result.retained_singular_values =
        information.retained_singular_values;

    const Eigen::Index levels =
        static_cast<Eigen::Index>(count);

    const Eigen::Index resolved =
        information.resolved_basis.cols();

    const Eigen::Index unresolved =
        information.unresolved_basis.cols();

    if (result.retained_singular_values.size() != resolved)
        throw std::invalid_argument(
            "Inconsistent information decomposition");

    result.resolved_sensitivity =
        Eigen::MatrixXd::Zero(levels, resolved);

    result.unresolved_sensitivity =
        Eigen::MatrixXd::Zero(levels, unresolved);

    result.maximum_nominal_mass_defect =
        recorded.summary.max_mass_defect;

    propagate_basis(
        problem,
        nominal,
        information.resolved_basis,
        configuration.parameters,
        result.resolved_sensitivity,
        result.maximum_tangent_residual,
        result.maximum_tangent_mass_defect);

    propagate_basis(
        problem,
        nominal,
        information.unresolved_basis,
        configuration.parameters,
        result.unresolved_sensitivity,
        result.maximum_tangent_residual,
        result.maximum_tangent_mass_defect);

    result.time.resize(count);
    result.nominal_outlet_pressure.resize(count);
    result.measurement_unit_radius.resize(count);
    result.unresolved_unit_radius.resize(count);
    result.unresolved_influence.resize(count);

    for (Eigen::Index k = 0; k < levels; ++k) {
        const auto index = static_cast<std::size_t>(k);

        result.time[index] =
            static_cast<double>(k) *
            configuration.parameters.dt;

        result.nominal_outlet_pressure[index] =
            recorded.summary.outlet_pressure[index];

        double measured_squared = 0.0;

        for (Eigen::Index j = 0; j < resolved; ++j) {
            const double singular =
                result.retained_singular_values[j];

            if (!(singular > 0.0) ||
                !std::isfinite(singular))
                throw std::runtime_error(
                    "Invalid retained singular value");

            const double contribution =
                result.resolved_sensitivity(k, j) / singular;

            measured_squared += contribution * contribution;
        }

        result.measurement_unit_radius[index] =
            std::sqrt(measured_squared);

        result.unresolved_unit_radius[index] =
            result.unresolved_sensitivity.row(k).norm();

        result.unresolved_influence[index] =
            result.unresolved_unit_radius[index] > 0.0;

        if (!std::isfinite(
                result.measurement_unit_radius[index]) ||
            !std::isfinite(
                result.unresolved_unit_radius[index]))
            throw std::runtime_error(
                "Nonfinite directional uncertainty");
    }

    return result;
}

} // namespace tpsf
