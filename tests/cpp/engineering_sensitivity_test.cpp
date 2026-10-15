#include "tpsf/measurement_uncertainty.hpp"
#include "tpsf/engineering_sensitivity.hpp"
#include "tpsf/dense_full_space_subproblem.hpp"
#include "tpsf/scheduled_forward.hpp"

#include <Eigen/Dense>

#include <algorithm>
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

void test_sensitivity(bool observed)
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
            "Incorrect unresolved dimension");

    const auto assembled =
        assemble_dense_full_space_subproblem(
            problem, solution.coordinates);

    const auto direction =
        quality.feasible_basis.col(0);

    require(
        (assembled.constraint_jacobian * direction).norm()
            < 1.0e-8,
        "Direction is not tangent to constraints");

    for (const auto quantity : {
        EngineeringQuantity::initial_linepack,
        EngineeringQuantity::terminal_linepack,
        EngineeringQuantity::minimum_pressure
    }) {
        const auto sensitivity =
            assess_engineering_sensitivity(
                problem,
                make_parameters(),
                solution.coordinates,
                quality,
                quantity);

        if (!sensitivity.differentiable) {
            require(
                quantity == EngineeringQuantity::minimum_pressure,
                "Unexpected nonsmooth engineering quantity");

            std::cout
                << "minimum pressure has "
                << sensitivity.active_minimum_count
                << " active locations; gradient not unique\n";

            continue;
        }

        constexpr double h = 2.0e-4;

        const double plus =
            evaluate_engineering_quantity(
                problem,
                make_parameters(),
                solution.coordinates + h*direction,
                quantity);

        const double minus =
            evaluate_engineering_quantity(
                problem,
                make_parameters(),
                solution.coordinates - h*direction,
                quantity);

        const double independent =
            (plus - minus)/(2.0*h);

        const double projected =
            sensitivity.gradient.dot(direction);

        const double tolerance =
            1.0e-3 *
            std::max({1.0,
                      std::abs(independent),
                      std::abs(projected)});

        require(
            std::abs(independent - projected) <= tolerance,
            "Directional sensitivity validation failed");

        require(
            std::abs(
                sensitivity.feasible_sensitivity_norm -
                std::abs(projected)) <= tolerance,
            "Incorrect feasible sensitivity projection");

        if (observed) {
            require(
                sensitivity.unresolved_sensitivity_norm <
                    1.0e-12,
                "Observed case has unresolved sensitivity");
        } else {
            require(
                std::abs(
                    sensitivity.unresolved_sensitivity_norm -
                    std::abs(projected)) <= tolerance,
                "Incorrect unresolved sensitivity projection");
        }

        std::cout
            << "observed=" << observed
            << " quantity=" << static_cast<int>(quantity)
            << " value=" << sensitivity.value
            << " feasible_norm="
            << sensitivity.feasible_sensitivity_norm
            << " unresolved_norm="
            << sensitivity.unresolved_sensitivity_norm
            << " derivative_check="
            << std::abs(independent-projected)
            << '\n';
    }
}


void test_measurement_uncertainty_integration(bool observed)
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
            "Expected one feasible state direction");

    require(quality.unresolved_directions ==
                static_cast<std::size_t>(observed ? 0 : 1),
            "Incorrect observability classification");

    const auto assembled =
        assemble_dense_full_space_subproblem(
            problem, solution.coordinates);

    const Eigen::MatrixXd S =
        assembled.measurement_jacobian *
        quality.feasible_basis;

    for (const auto quantity : {
        EngineeringQuantity::initial_linepack,
        EngineeringQuantity::terminal_linepack,
        EngineeringQuantity::minimum_pressure
    }) {
        const auto sensitivity =
            assess_engineering_sensitivity(
                problem,
                make_parameters(),
                solution.coordinates,
                quality,
                quantity);

        require(sensitivity.differentiable,
                "Integration fixture has nonsmooth quantity");

        const auto result =
            assess_engineering_uncertainty(
                problem,
                solution.coordinates,
                quality,
                sensitivity);

        const Eigen::VectorXd g =
            quality.feasible_basis.transpose() *
            sensitivity.gradient;

        const auto direct =
            propagate_measurement_uncertainty(S, g);

        require(result.status == direct.status,
                "Adapter and direct uncertainty disagree");

        require(result.information_rank ==
                    direct.information_rank,
                "Adapter and direct information ranks disagree");

        require(result.unresolved_directions ==
                    direct.unresolved_directions,
                "Adapter and direct nullspaces disagree");

        if (observed) {
            require(
                result.status == UncertaintyStatus::finite,
                "Observed quantity should have finite uncertainty");

            require(
                std::isfinite(result.standard_deviation),
                "Observed uncertainty is not finite");

            require(
                result.standard_deviation >= 0.0,
                "Negative standard deviation");

            const double tolerance =
                1.0e-10 *
                std::max(1.0, direct.standard_deviation);

            require(
                std::abs(
                    result.standard_deviation -
                    direct.standard_deviation) <= tolerance,
                "Adapter and direct uncertainty values disagree");
        } else {
            require(
                sensitivity.unresolved_sensitivity_norm > 1.0e-6,
                "Fixture quantity is not unresolved-sensitive");

            require(
                result.status == UncertaintyStatus::unbounded,
                "Unobserved quantity must be unbounded");

            require(
                !result.finite_data_only_bound,
                "Unobserved quantity incorrectly has finite bound");

            require(
                std::isinf(result.standard_deviation),
                "Unbounded uncertainty must be infinite");
        }

        std::cout
            << "E5.3 integration observed=" << observed
            << " quantity=" << static_cast<int>(quantity)
            << " information_rank=" << result.information_rank
            << " unresolved=" << result.unresolved_directions
            << " sigma=" << result.standard_deviation
            << '\n';
    }
}

} // namespace

int main()
{
    try {
        test_sensitivity(false);
        test_sensitivity(true);
        test_measurement_uncertainty_integration(false);
        test_measurement_uncertainty_integration(true);

        std::cout
            << "PASS E5.2 engineering quantity sensitivity\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E5.2: " << error.what() << '\n';

        return 1;
    }
}
