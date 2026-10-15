#include "tpsf/statefinding_quality.hpp"
#include "tpsf/dense_full_space_subproblem.hpp"
#include "tpsf/scheduled_forward.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

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

void test_observability()
{
    using namespace tpsf;

    for (const bool observed : {false, true}) {
        auto problem = make_problem(observed);

        DenseSQPResult solution;
        solution.coordinates =
            Eigen::VectorXd::Zero(
                static_cast<Eigen::Index>(
                    problem.free_variable_count()));

        solution.termination =
            DenseSQPTermination::converged;

        const auto quality =
            assess_statefinding_quality(problem, solution);

        require(quality.variables == 15,
                "Unexpected free-variable count");

        require(quality.feasible_dimension == 1,
                "Expected one feasible direction");

        require(quality.constraint_rank == 14,
                "Expected one independent constraint");

        require(quality.measurements ==
                    static_cast<std::size_t>(observed ? 1 : 0),
                "Incorrect measurement count");

        require(quality.measurement_rank ==
                    static_cast<std::size_t>(observed ? 1 : 0),
                "Incorrect measurement information rank");

        require(quality.unresolved_directions ==
                    static_cast<std::size_t>(observed ? 0 : 1),
                "Incorrect unresolved dimension");

        require(quality.locally_identifiable == observed,
                "Incorrect identifiability classification");

        require(quality.constraint_l2 < 1.0e-8,
                "Reference constraints not satisfied");

        require(quality.measurement_l2 < 1.0e-8,
                "Reference measurements not satisfied");

        const auto& Z = quality.feasible_basis;

        const auto assembled =
            assemble_dense_full_space_subproblem(
                problem, solution.coordinates);

        require(
            (assembled.constraint_jacobian * Z).norm() <
                1.0e-8,
            "Feasible basis does not satisfy A*Z=0");

        require(
            (Z.transpose()*Z -
             Eigen::MatrixXd::Identity(
                 Z.cols(), Z.cols())).norm() < 1.0e-10,
            "Feasible basis is not orthonormal");

        std::cout
            << "observed=" << observed
            << " feasible=" << quality.feasible_dimension
            << " measurement_rank=" << quality.measurement_rank
            << " unresolved=" << quality.unresolved_directions
            << '\n';
    }
}

} // namespace

int main()
{
    try {
        test_observability();
        std::cout
            << "PASS E5.1 statefinding quality and observability\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E5.1: " << error.what() << '\n';
        return 1;
    }
}
