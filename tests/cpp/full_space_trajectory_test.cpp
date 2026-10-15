#include "tpsf/full_space_trajectory.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_close(
    double actual,
    double expected,
    double tolerance,
    const std::string& message)
{
    require(
        std::isfinite(actual) &&
        std::isfinite(expected) &&
        std::abs(actual-expected) <= tolerance,
        message);
}

struct Fixture {
    static constexpr int n = 5;
    static constexpr int steps = 3;

    pipe_sim::TransientParameters parameters{
        100.0, 1.0, 0.65,
        0.5, 0.19634954084936207,
        288.15, 0.9, 518.3,
        1.1e-5, 1.0e-5
    };

    tpsf::FullSpaceTrajectory trajectory;

    Fixture()
    {
        trajectory.n = n;

        for (int k = 0; k <= steps; ++k) {
            std::vector<double> state(2*n+2);

            state[0] = 5.0e6 + 1000.0*k;
            state[1] = 10.0 + 0.1*k;

            for (int i = 0; i < n; ++i) {
                state[2*i+2] =
                    10.0 + 0.2*i + 0.1*k;

                state[2*i+3] =
                    5.0e6 - 20000.0*(i+1)
                    + 1000.0*k;
            }

            const tpsf::PhysicalBoundaryInputs inputs{
                state[0],
                11.0 + 0.15*k
            };

            const auto physical =
                tpsf::extract_physical_boundary_channels(
                    n, state, inputs);

            trajectory.states.push_back(state);
            trajectory.inputs.push_back(inputs);
            trajectory.reconstructed.push_back(physical);
        }
    }
};

void test_dimensions_and_ordering()
{
    const Fixture f;

    const auto result =
        tpsf::assemble_full_space_trajectory_constraints(
            f.trajectory, f.parameters);

    require(
        result.dynamics.size() == Fixture::steps,
        "Incorrect number of dynamics blocks");

    require(
        result.observation_links.size() == Fixture::steps+1,
        "Incorrect number of observation blocks");

    const auto combined = result.combined();

    const std::size_t expected_size =
        Fixture::steps*(2*Fixture::n+2)
        + 4*(Fixture::steps+1);

    require(
        combined.size() == expected_size,
        "Incorrect combined constraint dimension");

    std::size_t offset = 0;

    // Compare every dynamics block with the independently
    // assembled single-step residual.
    for (int k = 0; k < Fixture::steps; ++k) {
        const auto step =
            tpsf::assemble_full_space_step_constraints(
                Fixture::n,
                f.trajectory.states[k],
                f.trajectory.states[k+1],
                f.trajectory.inputs[k],
                f.trajectory.inputs[k+1],
                f.trajectory.reconstructed[k+1],
                f.parameters);

        require(
            step.dynamics.size() ==
                result.dynamics[k].size(),
            "Dynamics block dimension mismatch");

        for (std::size_t i = 0;
             i < step.dynamics.size(); ++i) {
            require_close(
                result.dynamics[k][i],
                step.dynamics[i],
                0.0,
                "Dynamics block differs from step residual");

            require_close(
                combined[offset++],
                step.dynamics[i],
                0.0,
                "Combined dynamics ordering incorrect");
        }
    }

    // Verify every observation-link block is placed after
    // all dynamics blocks, including the initial-time link.
    for (int k = 0; k <= Fixture::steps; ++k) {
        for (std::size_t c = 0; c < 4; ++c) {
            require_close(
                combined[offset++],
                result.observation_links[k][c],
                0.0,
                "Combined observation ordering incorrect");
        }
    }

    require(offset == combined.size(),
            "Combined constraint traversal incomplete");
}

void test_initial_link()
{
    Fixture f;

    const double perturbation = 1234.0;

    f.trajectory.reconstructed[0][0] += perturbation;

    const auto result =
        tpsf::assemble_full_space_trajectory_constraints(
            f.trajectory, f.parameters);

    require_close(
        result.observation_links[0][0],
        perturbation,
        1.0e-12,
        "Initial inlet-pressure link incorrect");

    for (const auto& dynamics : result.dynamics) {
        for (double value : dynamics) {
            require(std::isfinite(value),
                    "Nonfinite dynamics after link perturbation");
        }
    }

    const auto baseline =
        [&]() {
            Fixture original;
            return
                tpsf::assemble_full_space_trajectory_constraints(
                    original.trajectory,
                    original.parameters);
        }();

    for (std::size_t k = 0;
         k < result.dynamics.size(); ++k) {
        for (std::size_t i = 0;
             i < result.dynamics[k].size(); ++i) {
            require_close(
                result.dynamics[k][i],
                baseline.dynamics[k][i],
                0.0,
                "Reconstructed initial channel changed dynamics");
        }
    }

    for (std::size_t k = 0;
         k < result.observation_links.size(); ++k) {
        for (std::size_t c = 0; c < 4; ++c) {
            const double expected =
                (k == 0 && c == 0)
                    ? perturbation
                    : 0.0;

            require_close(
                result.observation_links[k][c],
                expected,
                1.0e-12,
                "Unexpected observation-link coupling");
        }
    }
}

void test_reconstructed_independence()
{
    Fixture f;

    const auto baseline =
        tpsf::assemble_full_space_trajectory_constraints(
            f.trajectory, f.parameters);

    const int k = 2;
    const int channel = 3;
    const double perturbation = 0.25;

    f.trajectory.reconstructed[k][channel] +=
        perturbation;

    const auto modified =
        tpsf::assemble_full_space_trajectory_constraints(
            f.trajectory, f.parameters);

    require(
        baseline.dynamics == modified.dynamics,
        "Reconstructed channel changed physical dynamics");

    for (std::size_t j = 0;
         j < modified.observation_links.size(); ++j) {
        for (std::size_t c = 0; c < 4; ++c) {
            const double expected =
                (j == static_cast<std::size_t>(k) &&
                 c == static_cast<std::size_t>(channel))
                    ? perturbation
                    : 0.0;

            require_close(
                modified.observation_links[j][c] -
                    baseline.observation_links[j][c],
                expected,
                1.0e-12,
                "Reconstructed-channel independence failed");
        }
    }
}

void test_input_independence()
{
    Fixture f;

    const auto baseline =
        tpsf::assemble_full_space_trajectory_constraints(
            f.trajectory, f.parameters);

    // Perturb the prescribed inlet pressure at time k=1.
    // This enters F_0 as its new inlet pressure and
    // affects C_1, but does not enter F_1 as an old
    // prescribed inlet pressure.
    constexpr int k = 1;
    constexpr double perturbation = 500.0;

    f.trajectory.inputs[k].inlet_pressure += perturbation;

    const auto modified =
        tpsf::assemble_full_space_trajectory_constraints(
            f.trajectory, f.parameters);

    require_close(
        modified.dynamics[0][0] -
            baseline.dynamics[0][0],
        -perturbation,
        0.0,
        "Prescribed inlet pressure not linked to F_0");

    for (std::size_t j = 1;
         j < modified.dynamics.size(); ++j) {
        require(
            modified.dynamics[j] == baseline.dynamics[j],
            "Prescribed inlet pressure affected wrong timestep");
    }

    require_close(
        modified.observation_links[k][0] -
            baseline.observation_links[k][0],
        -perturbation,
        0.0,
        "Prescribed inlet pressure not linked to C_1");
}

void test_invalid_dimensions()
{
    const Fixture f;

    auto expect_invalid = [&](tpsf::FullSpaceTrajectory trajectory) {
        bool rejected = false;

        try {
            (void)tpsf::assemble_full_space_trajectory_constraints(
                trajectory, f.parameters);
        }
        catch (const std::invalid_argument&) {
            rejected = true;
        }

        require(rejected,
                "Invalid trajectory dimensions accepted");
    };

    {
        auto t = f.trajectory;
        t.n = 1;
        expect_invalid(t);
    }

    {
        auto t = f.trajectory;
        t.states.pop_back();
        expect_invalid(t);
    }

    {
        auto t = f.trajectory;
        t.inputs.pop_back();
        expect_invalid(t);
    }

    {
        auto t = f.trajectory;
        t.reconstructed.pop_back();
        expect_invalid(t);
    }

    {
        auto t = f.trajectory;
        t.states[1].pop_back();
        expect_invalid(t);
    }
}

} // namespace

int main()
{
    try {
        test_dimensions_and_ordering();
        test_initial_link();
        test_reconstructed_independence();
        test_input_independence();
        test_invalid_dimensions();

        std::cout
            << "PASS full-space trajectory constraints\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
