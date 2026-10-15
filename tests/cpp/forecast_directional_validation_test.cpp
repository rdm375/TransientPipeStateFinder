#include "tpsf/forecast_directional_uncertainty.hpp"
#include "tpsf/measurement_information.hpp"
#include "tpsf/dense_full_space_subproblem.hpp"
#include "tpsf/forecast_continuation.hpp"
#include "tpsf/scheduled_forward.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace tpsf;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected,
          double atol, double rtol,
          const std::string& message)
{
    require(std::isfinite(actual) && std::isfinite(expected),
            message + ": nonfinite");

    require(std::abs(actual - expected) <=
                atol + rtol * std::max(
                    std::abs(actual), std::abs(expected)),
            message + ": mismatch");
}

std::vector<double> as_std(const Eigen::VectorXd& x)
{
    return {x.data(), x.data() + x.size()};
}

pipe_sim::TransientParameters parameters()
{
    return {
        100.0, 1.0, 0.65,
        0.5, 0.19634954084936207,
        288.15, 0.9, 518.3,
        1.1e-5, 1.0e-5
    };
}

Configuration forecast_configuration()
{
    return {
        2,
        10,
        parameters(),
        pipe_sim::DemandRamp{
            5.0e6, 10.1, 10.1, 10.0
        },
        pipe_sim::NewtonOptions{
            1.0e-10, 1.0e-12, 40, false
        }
    };
}

BoundarySchedule future_boundary()
{
    BoundarySchedule boundary;
    boundary.inlet_pressure.assign(11, 5.0e6);
    boundary.outlet_flow.assign(11, 10.1);
    return boundary;
}

FullSpaceProblem make_problem(bool observed)
{
    Configuration cfg{
        2,
        1,
        parameters(),
        pipe_sim::DemandRamp{
            5.0e6, 10.0, 10.1, 1.0
        },
        pipe_sim::NewtonOptions{
            1.0e-10, 1.0e-12, 40, false
        }
    };

    InitialState initial{
        {5.0e6, 5.0e6, 5.0e6},
        {10.0, 10.0},
        10.0
    };

    BoundarySchedule schedule{
        {5.0e6, 5.0e6},
        {10.0, 10.1}
    };

    const auto forward =
        run_forward_scheduled_recorded(
            cfg, initial, schedule);

    require(forward.summary.info == 0,
            "Fixture forward failed");

    FullSpaceTrajectory reference;
    reference.n = 2;

    for (const auto& snapshot :
         forward.trajectory.snapshots) {
        reference.states.push_back(snapshot.state);
        reference.inputs.push_back(snapshot.inputs);
        reference.reconstructed.push_back(
            snapshot.physical_channels);
    }

    FullSpaceProblemConfiguration problem_cfg;
    problem_cfg.reference = reference;
    problem_cfg.parameters = cfg.parameters;

    constexpr std::size_t width = 6;
    constexpr std::size_t block = width + 6;

    problem_cfg.variables.assign(
        2 * block,
        ProblemVariable{false, 1.0});

    for (std::size_t k = 0; k < 2; ++k) {
        const std::size_t offset = k * block;

        if (k > 0) {
            for (std::size_t j = 0; j < width; ++j) {
                const bool pressure =
                    j == 0 || (j >= 3 && j % 2 == 1);

                problem_cfg.variables[offset + j] = {
                    true,
                    pressure ? 1.0e6 : 10.0
                };
            }
        }

        for (std::size_t channel = 0;
             channel < 4; ++channel) {
            problem_cfg.variables[
                offset + width + 2 + channel] = {
                true,
                (channel == 0 || channel == 2)
                    ? 1.0e6 : 10.0
            };
        }
    }

    problem_cfg.variables[5] = {true, 1.0e6};

    const auto constraints =
        assemble_full_space_trajectory_constraints(
            reference, cfg.parameters);

    problem_cfg.constraint_scales.assign(
        constraints.combined().size(), 1.0);

    for (std::size_t k = 0; k < 2; ++k) {
        BoundaryMeasurement measurement;
        measurement.observed =
            reference.reconstructed[k];

        measurement.available = {
            false, false,
            observed && k == 1,
            false
        };

        measurement.sigma = {
            1.0e3, 0.1, 1.0e3, 0.1
        };

        problem_cfg.measurements.push_back(measurement);
    }

    return FullSpaceProblem(std::move(problem_cfg));
}

InitialState perturbed_initial(
    const InitialState& nominal,
    const Eigen::VectorXd& direction,
    double displacement)
{
    require(direction.size() == 6,
            "Incorrect physical direction width");

    InitialState result = nominal;

    result.pressure[0] += displacement * direction[0];
    result.inlet_flow += displacement * direction[1];

    for (int i = 0; i < 2; ++i) {
        result.flow[i] +=
            displacement * direction[2*i+2];

        result.pressure[i+1] +=
            displacement * direction[2*i+3];
    }

    return result;
}

void validate_direction(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& coordinates,
    const Eigen::VectorXd& basis_direction,
    const Eigen::VectorXd& predicted,
    const std::string& label)
{
    const auto configuration =
        forecast_configuration();

    const auto boundary =
        future_boundary();

    const auto reconstruction =
        problem.trajectory(as_std(coordinates));

    const auto nominal =
        initial_state_from_full_space(
            reconstruction,
            reconstruction.states.size() - 1);

    const auto mapped =
        problem.trajectory_direction(
            as_std(basis_direction));

    const auto& physical = mapped.states.back();

    const Eigen::VectorXd initial_direction =
        Eigen::Map<const Eigen::VectorXd>(
            physical.data(),
            static_cast<Eigen::Index>(physical.size()));

    require(std::abs(initial_direction[0]) <
                1.0e-8 * std::max(
                    1.0, initial_direction.norm()),
            "Direction changes prescribed inlet pressure");

    for (double h : {1.0e-2, 1.0e-3, 1.0e-4}) {
        const auto plus = perturbed_initial(
            nominal, initial_direction, h);

        const auto minus = perturbed_initial(
            nominal, initial_direction, -h);

        const auto fp = run_forward_scheduled(
            configuration, plus, boundary);

        const auto fm = run_forward_scheduled(
            configuration, minus, boundary);

        require(fp.info == 0 && fm.info == 0,
                "Perturbed nonlinear forecast failed");

        double max_error = 0.0;

        for (Eigen::Index k = 0;
             k < predicted.size(); ++k) {
            const auto index =
                static_cast<std::size_t>(k);

            const double derivative =
                (fp.outlet_pressure[index] -
                 fm.outlet_pressure[index]) / (2.0*h);

            max_error = std::max(
                max_error,
                std::abs(derivative - predicted[k]));

            near(
                derivative,
                predicted[k],
                2.0e-2,
                1.0e-3,
                label + ": nonlinear derivative");
        }

        std::cout << label
                  << " h=" << h
                  << " maximum derivative error="
                  << max_error << '\n';
    }
}

void test_directional_validation(bool observed)
{
    const auto problem = make_problem(observed);

    const Eigen::VectorXd coordinates =
        Eigen::VectorXd::Zero(
            static_cast<Eigen::Index>(
                problem.free_variable_count()));

    DenseSQPResult solution;
    solution.coordinates = coordinates;
    solution.termination =
        DenseSQPTermination::converged;

    const auto quality =
        assess_statefinding_quality(
            problem, solution);

    const auto assembled =
        assemble_dense_full_space_subproblem(
            problem, coordinates);

    const auto information =
        decompose_measurement_information(
            assembled.measurement_jacobian,
            quality.feasible_basis);

    const auto directional =
        assess_forecast_directional_uncertainty(
            problem,
            coordinates,
            information,
            forecast_configuration(),
            future_boundary());

    require(directional.time.size() == 11,
            "Incorrect forecast length");

    require(directional.maximum_tangent_residual < 1.0e-6,
            "Excessive tangent residual");

    require(directional.maximum_tangent_mass_defect < 1.0e-5,
            "Excessive tangent mass defect");

    require(information.resolved_basis.cols() ==
                (observed ? 1 : 0),
            "Unexpected resolved rank");

    require(information.unresolved_basis.cols() ==
                (observed ? 0 : 1),
            "Unexpected unresolved rank");

    for (Eigen::Index j = 0;
         j < information.resolved_basis.cols(); ++j) {
        validate_direction(
            problem,
            coordinates,
            information.resolved_basis.col(j),
            directional.resolved_sensitivity.col(j),
            "resolved");
    }

    for (Eigen::Index j = 0;
         j < information.unresolved_basis.cols(); ++j) {
        validate_direction(
            problem,
            coordinates,
            information.unresolved_basis.col(j),
            directional.unresolved_sensitivity.col(j),
            "unresolved");
    }

    for (Eigen::Index k = 0;
         k < directional.resolved_sensitivity.rows(); ++k) {
        const auto index =
            static_cast<std::size_t>(k);

        double weighted_squared = 0.0;
        double feasible_squared = 0.0;

        for (Eigen::Index j = 0;
             j < information.resolved_basis.cols(); ++j) {
            const double sensitivity =
                directional.resolved_sensitivity(k, j);

            const double singular =
                information.retained_singular_values[j];

            weighted_squared +=
                (sensitivity / singular) *
                (sensitivity / singular);

            feasible_squared +=
                sensitivity * sensitivity;
        }

        for (Eigen::Index j = 0;
             j < information.unresolved_basis.cols(); ++j) {
            const double sensitivity =
                directional.unresolved_sensitivity(k, j);

            feasible_squared +=
                sensitivity * sensitivity;
        }

        near(
            directional.measurement_unit_radius[index],
            std::sqrt(weighted_squared),
            1.0e-8, 1.0e-10,
            "Measurement radius");

        near(
            directional.unresolved_unit_radius[index],
            directional.unresolved_sensitivity.row(k).norm(),
            1.0e-8, 1.0e-10,
            "Unresolved radius");

        if (information.information_rank > 0) {
            const double smallest =
                information.retained_singular_values[
                    information.retained_singular_values.size()-1];

            const double old_bound =
                std::sqrt(feasible_squared) / smallest;

            require(
                directional.measurement_unit_radius[index]
                    <= old_bound + 1.0e-8,
                "New measurement radius exceeds old bound");
        }
    }

    std::cout
        << "PASS E5.5c.3 observed=" << observed
        << " resolved=" << information.information_rank
        << " unresolved=" << information.unresolved_dimension
        << '\n';
}

} // namespace

int main()
{
    try {
        test_directional_validation(false);
        test_directional_validation(true);

        std::cout
            << "PASS E5.5c.3 independent nonlinear validation\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E5.5c.3: "
            << error.what() << '\n';
        return 1;
    }
}
