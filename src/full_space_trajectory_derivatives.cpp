#include "tpsf/full_space_trajectory_derivatives.hpp"

#include "tpsf/full_space_tangent.hpp"
#include "tpsf/full_space_adjoint.hpp"
#include "tpsf/physical_boundary_channels.hpp"

#include <stdexcept>

namespace tpsf {
namespace {

void validate_trajectory(
    const FullSpaceTrajectory& trajectory)
{
    if (trajectory.n < 2 ||
        trajectory.states.size() < 2 ||
        trajectory.inputs.size() != trajectory.states.size() ||
        trajectory.reconstructed.size() != trajectory.states.size()) {
        throw std::invalid_argument(
            "Invalid full-space trajectory dimensions");
    }

    const auto nu =
        static_cast<std::size_t>(2*trajectory.n+2);

    for (const auto& state : trajectory.states) {
        if (state.size() != nu) {
            throw std::invalid_argument(
                "Invalid trajectory state dimension");
        }
    }
}

void validate_direction(
    const FullSpaceTrajectory& trajectory,
    const FullSpaceTrajectoryDirection& direction)
{
    if (direction.states.size() != trajectory.states.size() ||
        direction.inputs.size() != trajectory.inputs.size() ||
        direction.reconstructed.size() !=
            trajectory.reconstructed.size()) {
        throw std::invalid_argument(
            "Invalid trajectory direction dimensions");
    }

    const auto nu =
        static_cast<std::size_t>(2*trajectory.n+2);

    for (const auto& state : direction.states) {
        if (state.size() != nu) {
            throw std::invalid_argument(
                "Invalid trajectory tangent state dimension");
        }
    }
}

void validate_dual(
    const FullSpaceTrajectory& trajectory,
    const FullSpaceTrajectoryConstraints& dual)
{
    const auto steps = trajectory.steps();
    const auto nu =
        static_cast<std::size_t>(2*trajectory.n+2);

    if (dual.dynamics.size() != steps ||
        dual.observation_links.size() != steps+1) {
        throw std::invalid_argument(
            "Invalid trajectory dual block counts");
    }

    for (const auto& block : dual.dynamics) {
        if (block.size() != nu) {
            throw std::invalid_argument(
                "Invalid trajectory dynamics dual dimension");
        }
    }
}

void accumulate(
    std::vector<double>& destination,
    const std::vector<double>& contribution)
{
    if (destination.size() != contribution.size()) {
        throw std::invalid_argument(
            "Adjoint accumulation dimension mismatch");
    }

    for (std::size_t i = 0; i < destination.size(); ++i) {
        destination[i] += contribution[i];
    }
}

void accumulate(
    PhysicalBoundaryInputs& destination,
    const PhysicalBoundaryInputs& contribution)
{
    destination.inlet_pressure +=
        contribution.inlet_pressure;

    destination.outlet_flow +=
        contribution.outlet_flow;
}

} // namespace

FullSpaceTrajectoryConstraints
apply_full_space_trajectory_jvp(
    const FullSpaceTrajectory& trajectory,
    const FullSpaceTrajectoryDirection& direction,
    const pipe_sim::TransientParameters& parameters)
{
    validate_trajectory(trajectory);
    validate_direction(trajectory, direction);

    const auto steps = trajectory.steps();

    FullSpaceTrajectoryConstraints result;

    result.dynamics.resize(steps);
    result.observation_links.resize(steps+1);

    // The initial observation link is not covered by any
    // single-step derivative operator.
    const auto initial_channel_tangent =
        physical_boundary_channels_jvp(
            trajectory.n,
            direction.states[0],
            direction.inputs[0]);

    for (std::size_t c = 0; c < 4; ++c) {
        result.observation_links[0][c] =
            direction.reconstructed[0][c] -
            initial_channel_tangent[c];
    }

    // Each step contributes one dynamics block and the
    // observation link at its new time level.
    for (std::size_t k = 0; k < steps; ++k) {
        const auto step =
            apply_full_space_step_jvp(
                trajectory.n,
                trajectory.states[k],
                trajectory.states[k+1],
                direction.states[k],
                direction.states[k+1],
                direction.inputs[k],
                direction.inputs[k+1],
                direction.reconstructed[k+1],
                parameters);

        result.dynamics[k] = step.dynamics;

        result.observation_links[k+1] =
            step.observation_links;
    }

    return result;
}

FullSpaceTrajectoryAdjoint
apply_full_space_trajectory_vjp(
    const FullSpaceTrajectory& trajectory,
    const FullSpaceTrajectoryConstraints& dual,
    const pipe_sim::TransientParameters& parameters)
{
    validate_trajectory(trajectory);
    validate_dual(trajectory, dual);

    const auto count = trajectory.states.size();
    const auto nu =
        static_cast<std::size_t>(2*trajectory.n+2);

    FullSpaceTrajectoryAdjoint result;

    result.states.assign(
        count,
        std::vector<double>(nu, 0.0));

    result.inputs.resize(count);
    result.reconstructed.resize(count);

    // C_0 = b_0 - H(x_0, q_0).
    //
    // Handle the initial observation link explicitly.
    result.reconstructed[0] =
        dual.observation_links[0];

    BoundaryVector negative_initial{};

    for (std::size_t c = 0; c < 4; ++c) {
        negative_initial[c] =
            -dual.observation_links[0][c];
    }

    physical_boundary_channels_vjp(
        trajectory.n,
        negative_initial,
        result.states[0],
        result.inputs[0]);

    // Accumulate all shared-state and shared-input
    // contributions from adjacent timesteps.
    for (std::size_t k = 0; k < trajectory.steps(); ++k) {
        const auto step =
            apply_full_space_step_vjp(
                trajectory.n,
                trajectory.states[k],
                trajectory.states[k+1],
                dual.dynamics[k],
                dual.observation_links[k+1],
                parameters);

        accumulate(
            result.states[k],
            step.state_old);

        accumulate(
            result.states[k+1],
            step.state_new);

        accumulate(
            result.inputs[k],
            step.inputs_old);

        accumulate(
            result.inputs[k+1],
            step.inputs_new);

        result.reconstructed[k+1] =
            step.reconstructed_new;
    }

    return result;
}

} // namespace tpsf
