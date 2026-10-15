#pragma once

#include "tpsf/full_space_trajectory.hpp"
#include "tpsf/scheduled_forward.hpp"

#include <cstddef>

namespace tpsf {

// Convert a full-space state into the forward solver's
// initial-state representation.
InitialState initial_state_from_full_space(
    const FullSpaceTrajectory& trajectory,
    std::size_t time_level);

// Continue a reconstructed trajectory using the validated
// scheduled forward solver.
//
// The first future boundary values must match the
// reconstructed terminal boundary inputs.
RecordedForwardResult forecast_from_full_space(
    const Configuration& forecast_configuration,
    const FullSpaceTrajectory& reconstructed,
    const BoundarySchedule& future_boundary);

} // namespace tpsf
