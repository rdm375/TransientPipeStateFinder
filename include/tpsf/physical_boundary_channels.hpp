#pragma once

#include "tpsf/full_space_observations.hpp"

#include <span>

namespace tpsf {

struct PhysicalBoundaryInputs {
    double inlet_pressure = 0.0;
    double outlet_flow = 0.0;
};

BoundaryVector extract_physical_boundary_channels(
    int n,
    std::span<const double> state,
    const PhysicalBoundaryInputs& inputs);

BoundaryVector physical_boundary_channels_jvp(
    int n,
    std::span<const double> state_tangent,
    const PhysicalBoundaryInputs& input_tangent);

void physical_boundary_channels_vjp(
    int n,
    const BoundaryVector& channel_adjoint,
    std::span<double> state_adjoint,
    PhysicalBoundaryInputs& input_adjoint);

} // namespace tpsf
