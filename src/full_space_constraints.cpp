#include "tpsf/full_space_constraints.hpp"

#include <stdexcept>

namespace tpsf {

namespace {

void validate_state(int n, std::span<const double> state)
{
    if (n < 2 ||
        state.size() != static_cast<std::size_t>(2*n + 2)) {
        throw std::invalid_argument(
            "Invalid full-space state dimensions");
    }
}

} // namespace

std::vector<double> FullSpaceStepConstraints::combined() const
{
    std::vector<double> result;
    result.reserve(dynamics.size() + observation_links.size());

    result.insert(
        result.end(),
        dynamics.begin(),
        dynamics.end());

    result.insert(
        result.end(),
        observation_links.begin(),
        observation_links.end());

    return result;
}

FullSpaceStepConstraints assemble_full_space_step_constraints(
    int n,
    std::span<const double> state_old,
    std::span<const double> state_new,
    const PhysicalBoundaryInputs& inputs_old,
    const PhysicalBoundaryInputs& inputs_new,
    const BoundaryVector& reconstructed_new,
    const pipe_sim::TransientParameters& parameters)
{
    validate_state(n, state_old);
    validate_state(n, state_new);

    // Extract the old physical state in the exact ordering
    // required by the validated forward solver.
    std::vector<double> old_pressure(n + 1);
    std::vector<double> old_flow(n);

    old_pressure[0] = state_old[0];

    for (int i = 0; i < n; ++i) {
        old_flow[i] = state_old[2*i + 2];
        old_pressure[i + 1] = state_old[2*i + 3];
    }

    // Old inlet flow belongs to the previous physical state.
    // The outlet flow and inlet pressure are prescribed inputs.
    const pipe_sim::TransientBoundary boundary{
        .inlet_flow_old = state_old[1],
        .outlet_flow_old = inputs_old.outlet_flow,
        .outlet_flow_new = inputs_new.outlet_flow,
        .inlet_pressure_new = inputs_new.inlet_pressure
    };

    FullSpaceStepConstraints result;

    result.dynamics.resize(2*n + 2);

    // Use the validated forward solver's actual residual.
    pipe_sim::assemble_residual(
        n,
        state_new,
        old_pressure,
        old_flow,
        boundary,
        parameters,
        result.dynamics);

    // Extract physical boundary quantities.
    const BoundaryVector physical =
        extract_physical_boundary_channels(
            n,
            state_new,
            inputs_new);

    // Four independent observation-link constraints.
    for (std::size_t c = 0; c < 4; ++c) {
        result.observation_links[c] =
            reconstructed_new[c] - physical[c];
    }

    return result;
}

} // namespace tpsf
