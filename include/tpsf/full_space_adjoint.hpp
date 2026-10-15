#pragma once

#include "tpsf/full_space_constraints.hpp"

#include <span>
#include <vector>

namespace tpsf {

struct FullSpaceStepAdjoint {
    std::vector<double> state_old;
    std::vector<double> state_new;

    PhysicalBoundaryInputs inputs_old{};
    PhysicalBoundaryInputs inputs_new{};

    BoundaryVector reconstructed_new{};
};

// Exact transpose action for one full-space timestep.
//
// The dual variables correspond to:
//   dynamics          (2*n+2)
//   observation_links (4)
//
// Returns adjoints for all five independent variable groups.
//
// No finite differences or differentiation through Newton iterations.
FullSpaceStepAdjoint apply_full_space_step_vjp(
    int n,
    std::span<const double> state_old,
    std::span<const double> state_new,
    std::span<const double> dynamics_adjoint,
    const BoundaryVector& observation_adjoint,
    const pipe_sim::TransientParameters& parameters);

} // namespace tpsf
