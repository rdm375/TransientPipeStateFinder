#include "tpsf/full_space_trajectory.hpp"

#include <stdexcept>

namespace tpsf {

std::size_t FullSpaceTrajectory::steps() const
{
    if (states.empty())
        throw std::invalid_argument(
            "Trajectory must contain at least one state");

    return states.size()-1;
}

std::vector<double>
FullSpaceTrajectoryConstraints::combined() const
{
    std::size_t count = 4*observation_links.size();

    for (const auto& block : dynamics)
        count += block.size();

    std::vector<double> result;
    result.reserve(count);

    // All physical dynamics first.
    for (const auto& block : dynamics) {
        result.insert(
            result.end(), block.begin(), block.end());
    }

    // All observation links second.
    for (const auto& block : observation_links) {
        result.insert(
            result.end(), block.begin(), block.end());
    }

    return result;
}

FullSpaceTrajectoryConstraints
assemble_full_space_trajectory_constraints(
    const FullSpaceTrajectory& trajectory,
    const pipe_sim::TransientParameters& parameters)
{
    const int n = trajectory.n;

    if (n < 2)
        throw std::invalid_argument(
            "Trajectory requires at least two cells");

    const std::size_t count = trajectory.states.size();

    if (count < 2)
        throw std::invalid_argument(
            "Trajectory requires at least one timestep");

    if (trajectory.inputs.size() != count ||
        trajectory.reconstructed.size() != count) {
        throw std::invalid_argument(
            "Inconsistent trajectory time dimensions");
    }

    const std::size_t state_size =
        static_cast<std::size_t>(2*n+2);

    for (const auto& state : trajectory.states) {
        if (state.size() != state_size)
            throw std::invalid_argument(
                "Invalid trajectory state dimensions");
    }

    FullSpaceTrajectoryConstraints result;

    result.dynamics.reserve(count-1);
    result.observation_links.reserve(count);

    // --------------------------------------------------
    // Physical dynamics:
    // F_k(x_k, x_(k+1), q_k, q_(k+1)) = 0.
    // --------------------------------------------------

    for (std::size_t k = 0; k+1 < count; ++k) {
        const auto step =
            assemble_full_space_step_constraints(
                n,
                trajectory.states[k],
                trajectory.states[k+1],
                trajectory.inputs[k],
                trajectory.inputs[k+1],
                trajectory.reconstructed[k+1],
                parameters);

        result.dynamics.push_back(step.dynamics);
    }

    // --------------------------------------------------
    // Physical observation links at ALL time levels,
    // including k=0.
    //
    // C_k = b_k - H(x_k, q_k).
    // --------------------------------------------------

    for (std::size_t k = 0; k < count; ++k) {
        const auto physical =
            extract_physical_boundary_channels(
                n,
                trajectory.states[k],
                trajectory.inputs[k]);

        BoundaryVector link{};

        for (std::size_t c = 0; c < link.size(); ++c) {
            link[c] =
                trajectory.reconstructed[k][c]
                - physical[c];
        }

        result.observation_links.push_back(link);
    }

    return result;
}

} // namespace tpsf
