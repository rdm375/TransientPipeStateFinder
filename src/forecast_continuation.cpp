#include "tpsf/forecast_continuation.hpp"

#include <cmath>
#include <stdexcept>

namespace tpsf {

InitialState initial_state_from_full_space(
    const FullSpaceTrajectory& trajectory,
    std::size_t time_level)
{
    const int n = trajectory.n;

    if (n < 2)
        throw std::invalid_argument(
            "Forecast requires at least two cells");

    if (time_level >= trajectory.states.size())
        throw std::invalid_argument(
            "Invalid reconstruction time level");

    const auto& state = trajectory.states[time_level];

    if (state.size() != static_cast<std::size_t>(2*n+2))
        throw std::invalid_argument(
            "Invalid reconstructed state dimension");

    InitialState initial;

    initial.pressure.resize(n+1);
    initial.flow.resize(n);

    initial.pressure[0] = state[0];
    initial.inlet_flow = state[1];

    for (int i = 0; i < n; ++i) {
        initial.flow[i] = state[2*i+2];
        initial.pressure[i+1] = state[2*i+3];
    }

    for (double value : initial.pressure) {
        if (!std::isfinite(value) || value <= 0.0)
            throw std::invalid_argument(
                "Invalid reconstructed pressure");
    }

    if (!std::isfinite(initial.inlet_flow))
        throw std::invalid_argument(
            "Invalid reconstructed inlet flow");

    for (double value : initial.flow) {
        if (!std::isfinite(value))
            throw std::invalid_argument(
                "Invalid reconstructed internal flow");
    }

    return initial;
}

RecordedForwardResult forecast_from_full_space(
    const Configuration& forecast_configuration,
    const FullSpaceTrajectory& reconstructed,
    const BoundarySchedule& future_boundary)
{
    if (reconstructed.states.empty())
        throw std::invalid_argument(
            "Empty reconstructed trajectory");

    if (forecast_configuration.cells != reconstructed.n)
        throw std::invalid_argument(
            "Forecast spatial dimension mismatch");

    if (reconstructed.inputs.size() !=
        reconstructed.states.size())
        throw std::invalid_argument(
            "Reconstruction boundary history mismatch");

    if (future_boundary.inlet_pressure.empty() ||
        future_boundary.outlet_flow.empty())
        throw std::invalid_argument(
            "Empty forecast boundary schedule");

    const auto& terminal_inputs =
        reconstructed.inputs.back();

    // PhysicalBoundaryInputs ordering:
    // [inlet pressure, outlet mass flow].
    if (future_boundary.inlet_pressure.front() !=
        terminal_inputs.inlet_pressure)
        throw std::invalid_argument(
            "Forecast inlet pressure discontinuity");

    if (future_boundary.outlet_flow.front() !=
        terminal_inputs.outlet_flow)
        throw std::invalid_argument(
            "Forecast outlet flow discontinuity");

    const auto initial =
        initial_state_from_full_space(
            reconstructed,
            reconstructed.states.size()-1);

    return run_forward_scheduled_recorded(
        forecast_configuration,
        initial,
        future_boundary);
}

} // namespace tpsf
