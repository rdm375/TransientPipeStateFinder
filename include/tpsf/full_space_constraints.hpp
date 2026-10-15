#pragma once

#include "tpsf/full_space_observations.hpp"
#include "tpsf/physical_boundary_channels.hpp"

#include "pipe_sim/transient.hpp"

#include <span>
#include <vector>

namespace tpsf {

// Constraints for one physical timestep.
//
// Dynamics:
//   F(x_new, x_old, q_new, q_old) = 0
//
// Observation links:
//   C = b_new - H(x_new, q_new) = 0
//
// The four reconstructed channels b_new remain independent
// of the physical state and prescribed boundary inputs.

struct FullSpaceStepConstraints {
    std::vector<double> dynamics;
    BoundaryVector observation_links{};

    // Combined residual ordering:
    // [dynamics (2*n+2), observation links (4)]
    std::vector<double> combined() const;
};

FullSpaceStepConstraints assemble_full_space_step_constraints(
    int n,
    std::span<const double> state_old,
    std::span<const double> state_new,
    const PhysicalBoundaryInputs& inputs_old,
    const PhysicalBoundaryInputs& inputs_new,
    const BoundaryVector& reconstructed_new,
    const pipe_sim::TransientParameters& parameters);

} // namespace tpsf
