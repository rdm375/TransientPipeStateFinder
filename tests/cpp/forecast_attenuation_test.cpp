#include "tpsf/forecast_attenuation.hpp"
#include "tpsf/pipeline_decision.hpp"
#include "tpsf/scheduled_forward.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

tpsf::FullSpaceProblem make_problem(bool observed)
{
    using namespace tpsf;

    Configuration forward_cfg{
        2,
        1,
        pipe_sim::TransientParameters{
            100.0,
            1.0,
            0.65,
            0.5,
            0.19634954084936207,
            288.15,
            0.9,
            518.3,
            1.1e-5,
            1.0e-5
        },
        pipe_sim::DemandRamp{
            5.0e6, 10.0, 10.1, 1.0
        },
        pipe_sim::NewtonOptions{
            1.0e-8, 1.0e-10, 30, false
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
            forward_cfg, initial, schedule);

    require(forward.summary.info == 0,
            "Forward integration failed");

    require(forward.summary.completed_steps == 1,
            "Incomplete forward trajectory");

    FullSpaceTrajectory reference;
    reference.n = 2;

    for (const auto& snapshot : forward.trajectory.snapshots) {
        reference.states.push_back(snapshot.state);
        reference.inputs.push_back(snapshot.inputs);
        reference.reconstructed.push_back(
            snapshot.physical_channels);
    }

    require(reference.states.size() == 2,
            "Expected two trajectory time levels");

    FullSpaceProblemConfiguration cfg;
    cfg.reference = reference;
    cfg.parameters = forward_cfg.parameters;

    constexpr std::size_t state_size = 6;
    constexpr std::size_t block_size = state_size + 6;

    cfg.variables.assign(
        2 * block_size,
        ProblemVariable{false, 1.0});

    for (std::size_t k = 0; k < 2; ++k) {
        const std::size_t offset = k * block_size;

        // Evolving physical state is free after time zero.
        if (k > 0) {
            for (std::size_t j = 0; j < state_size; ++j) {
                const bool pressure =
                    j == 0 || (j >= 3 && j % 2 == 1);

                cfg.variables[offset + j] = {
                    true,
                    pressure ? 1.0e6 : 10.0
                };
            }
        }

        // Reconstructed boundary channels are free.
        for (std::size_t channel = 0; channel < 4; ++channel) {
            cfg.variables[offset + state_size + 2 + channel] = {
                true,
                (channel == 0 || channel == 2)
                    ? 1.0e6 : 10.0
            };
        }
    }

    // One additional independent initial-state variable.
    cfg.variables[5] = {true, 1.0e6};

    const auto constraints =
        assemble_full_space_trajectory_constraints(
            reference, cfg.parameters);

    cfg.constraint_scales.assign(
        constraints.combined().size(), 1.0);

    for (std::size_t k = 0; k < 2; ++k) {
        BoundaryMeasurement measurement;

        measurement.observed = reference.reconstructed[k];

        measurement.available = {
            false,
            false,
            observed && k == 1,
            false
        };

        measurement.sigma = {
            1.0e3, 0.1, 1.0e3, 0.1
        };

        cfg.measurements.push_back(measurement);
    }

    return FullSpaceProblem(std::move(cfg));
}



pipe_sim::TransientParameters make_parameters()
{
    return {
        100.0, 1.0, 0.65,
        0.5, 0.19634954084936207,
        288.15, 0.9, 518.3,
        1.1e-5, 1.0e-5
    };
}

void test_tangent_forecast(bool observed)
{
    using namespace tpsf;

    const auto problem = make_problem(observed);

    DenseSQPResult solution;
    solution.coordinates = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(
            problem.free_variable_count()));

    solution.termination = DenseSQPTermination::converged;

    const auto quality =
        assess_statefinding_quality(problem, solution);

    require(quality.feasible_dimension == 1,
            "Expected one feasible direction");

    require(quality.unresolved_directions ==
                static_cast<std::size_t>(observed ? 0 : 1),
            "Unexpected unresolved dimension");

    Configuration forecast{
        2,
        10,
        make_parameters(),
        pipe_sim::DemandRamp{
            5.0e6, 10.1, 10.1, 10.0
        },
        pipe_sim::NewtonOptions{
            1.0e-8, 1.0e-10, 30, false
        }
    };

    BoundarySchedule future;
    future.inlet_pressure.assign(11, 5.0e6);
    future.outlet_flow.assign(11, 10.1);

    const auto result =
        assess_forecast_attenuation(
            problem, solution.coordinates,
            quality, forecast, future);

    require(result.points.size() == 11,
            "Incorrect forecast length");

    require(result.maximum_mass_defect < 1.0e-5,
            "Excessive mass defect");

    require(result.maximum_tangent_residual < 1.0e-6,
            "Excessive tangent equation residual");

    for (const auto& point : result.points) {
        require(
            std::isfinite(point.linepack) &&
            std::isfinite(point.outlet_pressure) &&
            std::isfinite(point.linepack_feasible) &&
            std::isfinite(point.outlet_pressure_feasible),
            "Nonfinite forecast output");

        if (observed) {
            require(
                point.linepack_unresolved == 0.0 &&
                point.outlet_pressure_unresolved == 0.0,
                "Unexpected unresolved forecast sensitivity");
        } else {
            require(
                std::abs(point.linepack_unresolved -
                         point.linepack_feasible) < 1.0e-8,
                "Incorrect unresolved linepack sensitivity");

            require(
                std::abs(point.outlet_pressure_unresolved -
                         point.outlet_pressure_feasible) < 1.0e-6,
                "Incorrect unresolved pressure sensitivity");
        }

        std::cout
            << "observed=" << observed
            << " t=" << point.time
            << " linepack=" << point.linepack
            << " outlet_pressure=" << point.outlet_pressure
            << " linepack_unresolved="
            << point.linepack_unresolved
            << " pressure_unresolved="
            << point.outlet_pressure_unresolved
            << '\n';
    }

    std::cout
        << "maximum tangent residual: "
        << result.maximum_tangent_residual << '\n';
}


void test_pipeline_decision(bool observed)
{
    using namespace tpsf;

    const auto problem = make_problem(observed);

    DenseSQPResult solution;
    solution.coordinates = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(
            problem.free_variable_count()));
    solution.termination = DenseSQPTermination::converged;

    const auto quality =
        assess_statefinding_quality(problem, solution);

    Configuration forecast{
        2,
        10,
        make_parameters(),
        pipe_sim::DemandRamp{
            5.0e6, 10.1, 10.1, 10.0
        },
        pipe_sim::NewtonOptions{
            1.0e-8, 1.0e-10, 30, false
        }
    };

    BoundarySchedule future;
    future.inlet_pressure.assign(11, 5.0e6);
    future.outlet_flow.assign(11, 10.1);

    PipelineDecisionOptions options;
    options.minimum_outlet_pressure_pa = 4.0e6;
    options.measurement_radius = 1.0;

    // Unresolved directions are initially unbounded.
    const auto result = assess_pipeline_decision(
        problem,
        solution.coordinates,
        quality,
        forecast,
        future,
        options);

    require(result.decision.points.size() == 11,
            "Incorrect decision forecast length");

    require(result.attenuation.maximum_tangent_residual
                < 1.0e-6,
            "Excessive decision tangent residual");

    if (!observed) {
        require(
            result.decision.points.front().classification ==
                DecisionClassification::indeterminate,
            "Unobserved state must not receive a false certificate");

        require(
            !result.decision.points.front().finite_bounds,
            "Unresolved uncertainty must remain unbounded");

        // Explicitly bound unresolved scaled coordinates.
        options.bound_unresolved = true;
        options.unresolved_radius = 0.01;

        const auto bounded = assess_pipeline_decision(
            problem,
            solution.coordinates,
            quality,
            forecast,
            future,
            options);

        require(bounded.decision.points.front().finite_bounds,
                "Explicit uncertainty bound was not applied");

        require(
            bounded.decision.points.front().uncertainty_radius
                > 0.0,
            "Expected positive forecast uncertainty");

    } else {
        require(quality.unresolved_directions == 0,
                "Expected fully observed feasible subspace");

        require(
            result.decision.classification ==
                DecisionClassification::certified_satisfactory,
            "Expected satisfactory observed forecast");

        for (const auto& point : result.decision.points)
            require(point.finite_bounds,
                    "Observed forecast must have finite bounds");
    }

    std::cout
        << "PASS E5.5b pipeline decision observed="
        << observed << '\n';
}

} // namespace

int main()
{
    try {
        test_tangent_forecast(false);
        test_tangent_forecast(true);
        test_pipeline_decision(false);
        test_pipeline_decision(true);

        std::cout
            << "PASS E5.4 exact tangent-linear forecast\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E5.4: " << error.what() << '\n';
        return 1;
    }
}
