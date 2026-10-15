#pragma once

#include "tpsf/full_space_constraints.hpp"

namespace tpsf {

// Exact directional derivative of one full-space timestep.
//
// The base point contains the old and new physical states.
// Tangent arguments represent independent variations of:
//   x_old, x_new, q_old, q_new, b_new.
//
// Returns:
//   [DF (2*n+2), DC (4)]
//
// No finite differences are used.
FullSpaceStepConstraints apply_full_space_step_jvp(
    int n,
    std::span<const double> state_old,
    std::span<const double> state_new,
    std::span<const double> tangent_old,
    std::span<const double> tangent_new,
    const PhysicalBoundaryInputs& tangent_inputs_old,
    const PhysicalBoundaryInputs& tangent_inputs_new,
    const BoundaryVector& tangent_reconstructed_new,
    const pipe_sim::TransientParameters& parameters);

} // namespace tpsf
