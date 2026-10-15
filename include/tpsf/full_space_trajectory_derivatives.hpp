#pragma once

#include "tpsf/full_space_trajectory.hpp"

namespace tpsf {

// Independent variations at every trajectory time level.
struct FullSpaceTrajectoryDirection {
    std::vector<std::vector<double>> states;
    std::vector<PhysicalBoundaryInputs> inputs;
    std::vector<BoundaryVector> reconstructed;
};

// Accumulated transpose contributions for every variable.
struct FullSpaceTrajectoryAdjoint {
    std::vector<std::vector<double>> states;
    std::vector<PhysicalBoundaryInputs> inputs;
    std::vector<BoundaryVector> reconstructed;
};

// Exact matrix-free trajectory tangent action.
FullSpaceTrajectoryConstraints
apply_full_space_trajectory_jvp(
    const FullSpaceTrajectory& trajectory,
    const FullSpaceTrajectoryDirection& direction,
    const pipe_sim::TransientParameters& parameters);

// Exact matrix-free trajectory transpose action.
FullSpaceTrajectoryAdjoint
apply_full_space_trajectory_vjp(
    const FullSpaceTrajectory& trajectory,
    const FullSpaceTrajectoryConstraints& dual,
    const pipe_sim::TransientParameters& parameters);

} // namespace tpsf
