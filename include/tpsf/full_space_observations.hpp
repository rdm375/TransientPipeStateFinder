#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace tpsf {

constexpr std::size_t boundary_channels = 4;

enum class BoundaryChannel : std::size_t {
    inlet_pressure = 0,
    inlet_flow = 1,
    outlet_pressure = 2,
    outlet_flow = 3
};

using BoundaryVector = std::array<double, boundary_channels>;

// These are independent variational objects.
struct FullSpaceObservation {
    BoundaryVector reconstructed{};
    BoundaryVector measured{};
    BoundaryVector scales{};
    BoundaryVector state_boundary{};
    bool included{true};
};

struct FullSpaceObservationResidual {
    BoundaryVector measurement{};
    BoundaryVector constraint{};
};

// Residuals:
// measurement[c] = (reconstructed[c] - measured[c]) / scales[c]
// constraint[c]  = reconstructed[c] - state_boundary[c]
FullSpaceObservationResidual evaluate_full_space_observation(
    const FullSpaceObservation& observation);

// Jacobian action with respect to reconstructed boundary variables
// and the boundary values extracted from the pipe state.
FullSpaceObservationResidual full_space_observation_jvp(
    const FullSpaceObservation& observation,
    const BoundaryVector& delta_reconstructed,
    const BoundaryVector& delta_state_boundary);

// Transpose-Jacobian action.
// Returns gradients with respect to reconstructed boundary variables
// and state boundary values, respectively.
struct FullSpaceObservationAdjoint {
    BoundaryVector reconstructed{};
    BoundaryVector state_boundary{};
};

FullSpaceObservationAdjoint full_space_observation_vjp(
    const FullSpaceObservation& observation,
    const FullSpaceObservationResidual& dual);

// Measurement objective, excluding physical constraints.
double full_space_measurement_cost(
    const FullSpaceObservation& observation);

} // namespace tpsf
