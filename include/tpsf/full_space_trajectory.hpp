#pragma once

#include "tpsf/full_space_constraints.hpp"

#include <cstddef>
#include <vector>

namespace tpsf {

// Independent optimization variables at every time level.
//
// Time grid:
//   t_k = k * parameters.dt, k = 0,...,N.
struct FullSpaceTrajectory {
    int n = 0;

    // Physical state at each time level.
    // states[k].size() == 2*n+2
    std::vector<std::vector<double>> states;

    // Independently prescribed physical boundary inputs.
    std::vector<PhysicalBoundaryInputs> inputs;

    // Independently reconstructed four-channel boundaries.
    std::vector<BoundaryVector> reconstructed;

    std::size_t steps() const;
};

// Residual ordering:
//
// [ F_0, ..., F_(N-1), C_0, ..., C_N ]
//
// F_k contains 2*n+2 equations.
// C_k contains four boundary-channel equations.
struct FullSpaceTrajectoryConstraints {
    std::vector<std::vector<double>> dynamics;
    std::vector<BoundaryVector> observation_links;

    std::vector<double> combined() const;
};

FullSpaceTrajectoryConstraints
assemble_full_space_trajectory_constraints(
    const FullSpaceTrajectory& trajectory,
    const pipe_sim::TransientParameters& parameters);

} // namespace tpsf
