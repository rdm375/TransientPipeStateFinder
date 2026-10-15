#include "tpsf/physical_boundary_channels.hpp"

#include <stdexcept>

namespace tpsf {

namespace {

void validate_state(int n, std::size_t size)
{
    if (n < 2 ||
        size != static_cast<std::size_t>(2*n + 2)) {
        throw std::invalid_argument(
            "Invalid physical boundary state dimensions");
    }
}

} // namespace

BoundaryVector extract_physical_boundary_channels(
    int n,
    std::span<const double> state,
    const PhysicalBoundaryInputs& inputs)
{
    validate_state(n, state.size());

    return {
        inputs.inlet_pressure,
        state[1],
        state[2*n + 1],
        inputs.outlet_flow
    };
}

BoundaryVector physical_boundary_channels_jvp(
    int n,
    std::span<const double> state_tangent,
    const PhysicalBoundaryInputs& input_tangent)
{
    validate_state(n, state_tangent.size());

    return {
        input_tangent.inlet_pressure,
        state_tangent[1],
        state_tangent[2*n + 1],
        input_tangent.outlet_flow
    };
}

void physical_boundary_channels_vjp(
    int n,
    const BoundaryVector& channel_adjoint,
    std::span<double> state_adjoint,
    PhysicalBoundaryInputs& input_adjoint)
{
    validate_state(n, state_adjoint.size());

    state_adjoint[1] += channel_adjoint[1];
    state_adjoint[2*n + 1] += channel_adjoint[2];

    input_adjoint.inlet_pressure += channel_adjoint[0];
    input_adjoint.outlet_flow += channel_adjoint[3];
}

} // namespace tpsf
