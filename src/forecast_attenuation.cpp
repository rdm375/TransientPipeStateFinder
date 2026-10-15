#include "tpsf/forecast_attenuation.hpp"

#include "tpsf/full_space_tangent.hpp"
#include "pipe_sim/physics.hpp"

#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace tpsf {
namespace {

std::vector<double> as_std(const Eigen::VectorXd& u)
{
    return std::vector<double>(u.data(), u.data()+u.size());
}

Eigen::VectorXd as_eigen(const std::vector<double>& x)
{
    return Eigen::Map<const Eigen::VectorXd>(
        x.data(), static_cast<Eigen::Index>(x.size()));
}

std::vector<double> pressures(
    const std::vector<double>& state, int n)
{
    std::vector<double> result(n+1);
    result[0] = state[0];

    for (int i = 0; i < n; ++i)
        result[i+1] = state[2*i+3];

    return result;
}

double linepack_direction(
    const std::vector<double>& tangent,
    int n,
    const pipe_sim::TransientParameters& par)
{
    // Exact derivative of trapezoidal linepack for the
    // constant-z isothermal EOS used by the forward solver.
    const double crho =
        pipe_sim::density_derivative_cz(
            par.temperature, par.z, par.gas_constant);

    const auto dp = pressures(tangent, n);

    double weighted = 0.5*(dp.front()+dp.back());

    for (int i = 1; i < n; ++i)
        weighted += dp[i];

    return par.area*par.dx*crho*weighted;
}

ForecastTangentHistory propagate_internal(
    const FullSpaceTrajectory& nominal,
    const Eigen::VectorXd& initial_direction,
    const pipe_sim::TransientParameters& par)
{
    const int n = nominal.n;
    const Eigen::Index width = 2*n+2;

    if (initial_direction.size() != width)
        throw std::invalid_argument(
            "Incorrect initial tangent dimension");

    ForecastTangentHistory result;
    result.states.push_back(initial_direction);

    const PhysicalBoundaryInputs zero_inputs{};
    const BoundaryVector zero_channels{};

    const std::vector<double> zero(
        static_cast<std::size_t>(width), 0.0);

    for (std::size_t k = 0; k < nominal.steps(); ++k) {
        const auto& old_state = nominal.states[k];
        const auto& new_state = nominal.states[k+1];

        const auto old_tangent =
            as_std(result.states.back());

        // r = B_k * delta_x_old
        const auto rhs_action = apply_full_space_step_jvp(
            n,
            old_state,
            new_state,
            old_tangent,
            zero,
            zero_inputs,
            zero_inputs,
            zero_channels,
            par);

        Eigen::MatrixXd A(width, width);

        // Construct A_k from exact JVP actions.
        // This is a dense validation implementation.
        for (Eigen::Index j = 0; j < width; ++j) {
            auto unit = zero;
            unit[static_cast<std::size_t>(j)] = 1.0;

            const auto column =
                apply_full_space_step_jvp(
                    n,
                    old_state,
                    new_state,
                    zero,
                    unit,
                    zero_inputs,
                    zero_inputs,
                    zero_channels,
                    par);

            A.col(j) = as_eigen(column.dynamics);
        }

        const Eigen::VectorXd rhs =
            -as_eigen(rhs_action.dynamics);

        const Eigen::FullPivLU<Eigen::MatrixXd> lu(A);

        if (!lu.isInvertible())
            throw std::runtime_error(
                "Singular implicit forecast tangent Jacobian");

        const Eigen::VectorXd next = lu.solve(rhs);

        if (!next.allFinite())
            throw std::runtime_error(
                "Nonfinite forecast tangent state");

        const double residual =
            (A*next-rhs).lpNorm<Eigen::Infinity>();

        result.maximum_equation_residual =
            std::max(result.maximum_equation_residual, residual);

        result.states.push_back(next);
    }

    return result;
}

void accumulate_basis(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& coordinates,
    const FullSpaceTrajectory& reconstruction,
    const FullSpaceTrajectory& forecast,
    const Eigen::MatrixXd& basis,
    const pipe_sim::TransientParameters& par,
    std::vector<double>& linepack_squared,
    std::vector<double>& pressure_squared,
    double& maximum_residual)
{
    const Eigen::Index width = 2*forecast.n+2;

    for (Eigen::Index j = 0; j < basis.cols(); ++j) {
        // FullSpaceProblem uses an affine map from scaled
        // coordinates to physical trajectory variables.
        // Its exact directional derivative is therefore
        // trajectory(u+v)-trajectory(u), with no truncation
        // error or numerical differentiation of the model.
        const auto direction =
            problem.trajectory_direction(
                as_std(basis.col(j)));

        const Eigen::VectorXd initial =
            as_eigen(direction.states.back());

        // The prescribed inlet pressure is fixed during
        // forecast differentiation.
        const double inlet_tolerance =
            1.0e-10*std::max(1.0, initial.norm());

        if (std::abs(initial[0]) > inlet_tolerance)
            throw std::invalid_argument(
                "Tangent changes prescribed forecast inlet pressure");

        const auto history =
            propagate_forecast_tangent(
                forecast, initial, par);

        maximum_residual =
            std::max(maximum_residual,
                     history.maximum_equation_residual);

        for (std::size_t k = 0;
             k < history.states.size(); ++k) {
            const auto tangent =
                as_std(history.states[k]);

            const double d_mass =
                linepack_direction(tangent, forecast.n, par);

            const double d_outlet =
                tangent[static_cast<std::size_t>(
                    2*forecast.n+1)];

            linepack_squared[k] += d_mass*d_mass;
            pressure_squared[k] += d_outlet*d_outlet;
        }
    }
}

} // namespace

ForecastTangentHistory propagate_forecast_tangent(
    const FullSpaceTrajectory& nominal_forecast,
    const Eigen::VectorXd& initial_direction,
    const pipe_sim::TransientParameters& parameters)
{
    auto result = propagate_internal(
        nominal_forecast, initial_direction, parameters);

    const double dt = parameters.dt;
    const double theta = parameters.theta;

    if (!std::isfinite(dt) || dt <= 0.0 ||
        !std::isfinite(theta) ||
        theta < 0.0 || theta > 1.0)
        throw std::invalid_argument(
            "Invalid tangent conservation parameters");

    const int n = nominal_forecast.n;

    for (std::size_t k = 0;
         k + 1 < result.states.size(); ++k) {

        const auto old_tangent = as_std(result.states[k]);
        const auto new_tangent = as_std(result.states[k+1]);

        const double old_mass =
            linepack_direction(old_tangent, n, parameters);

        const double new_mass =
            linepack_direction(new_tangent, n, parameters);

        // State index 1 is inlet mass flow.
        // Future outlet-flow perturbations are zero.
        const double old_inlet = old_tangent[1];
        const double new_inlet = new_tangent[1];

        const double defect =
            (new_mass-old_mass)/dt
            - theta*new_inlet
            - (1.0-theta)*old_inlet;

        result.maximum_mass_balance_defect =
            std::max(
                result.maximum_mass_balance_defect,
                std::abs(defect));
    }

    return result;
}

ForecastAttenuationResult assess_forecast_attenuation(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& coordinates,
    const StatefindingQuality& quality,
    const Configuration& configuration,
    const BoundarySchedule& boundary)
{
    const Eigen::Index n =
        static_cast<Eigen::Index>(
            problem.free_variable_count());

    if (coordinates.size() != n ||
        !coordinates.allFinite() ||
        quality.feasible_basis.rows() != n ||
        quality.unresolved_basis.rows() != n)
        throw std::invalid_argument(
            "Invalid forecast attenuation dimensions");

    if (configuration.steps < 1)
        throw std::invalid_argument(
            "Forecast must contain at least one time step");

    const std::size_t count =
        static_cast<std::size_t>(configuration.steps+1);

    if (boundary.inlet_pressure.size() != count ||
        boundary.outlet_flow.size() != count)
        throw std::invalid_argument(
            "Incorrect forecast boundary schedule length");

    const auto reconstruction =
        problem.trajectory(as_std(coordinates));

    const auto recorded =
        forecast_from_full_space(
            configuration, reconstruction, boundary);

    if (recorded.summary.info != 0 ||
        recorded.summary.completed_steps != configuration.steps)
        throw std::runtime_error(
            "Nominal forecast failed");

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
        recorded.summary.linepack.size() != count ||
        recorded.summary.outlet_pressure.size() != count)
        throw std::runtime_error(
            "Incomplete nominal forecast recording");

    ForecastAttenuationResult result;
    result.points.resize(count);

    result.maximum_mass_defect =
        recorded.summary.max_mass_defect;

    std::vector<double> fl(count, 0.0);
    std::vector<double> fp(count, 0.0);
    std::vector<double> ul(count, 0.0);
    std::vector<double> up(count, 0.0);

    accumulate_basis(
        problem, coordinates, reconstruction, nominal,
        quality.feasible_basis,
        configuration.parameters,
        fl, fp, result.maximum_tangent_residual);

    accumulate_basis(
        problem, coordinates, reconstruction, nominal,
        quality.unresolved_basis,
        configuration.parameters,
        ul, up, result.maximum_tangent_residual);

    for (std::size_t k = 0; k < count; ++k) {
        auto& point = result.points[k];

        point.time =
            static_cast<double>(k) *
            configuration.parameters.dt;

        point.linepack = recorded.summary.linepack[k];
        point.outlet_pressure =
            recorded.summary.outlet_pressure[k];

        point.linepack_feasible = std::sqrt(fl[k]);
        point.outlet_pressure_feasible = std::sqrt(fp[k]);
        point.linepack_unresolved = std::sqrt(ul[k]);
        point.outlet_pressure_unresolved = std::sqrt(up[k]);
    }

    return result;
}

} // namespace tpsf
