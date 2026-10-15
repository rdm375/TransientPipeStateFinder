#include "tpsf/full_space_observations.hpp"

#include <cmath>
#include <stdexcept>

namespace tpsf {

namespace {

void validate(const FullSpaceObservation& observation)
{
    for (std::size_t c = 0; c < boundary_channels; ++c) {
        if (!std::isfinite(observation.reconstructed[c]) ||
            !std::isfinite(observation.measured[c]) ||
            !std::isfinite(observation.state_boundary[c]) ||
            !std::isfinite(observation.scales[c]) ||
            observation.scales[c] <= 0.0)
            throw std::invalid_argument(
                "Invalid full-space observation");
    }
}

} // namespace

FullSpaceObservationResidual evaluate_full_space_observation(
    const FullSpaceObservation& observation)
{
    validate(observation);

    FullSpaceObservationResidual result;

    for (std::size_t c = 0; c < boundary_channels; ++c) {
        result.constraint[c] =
            observation.reconstructed[c] -
            observation.state_boundary[c];

        if (observation.included) {
            result.measurement[c] =
                (observation.reconstructed[c] -
                 observation.measured[c]) /
                observation.scales[c];
        }
    }

    return result;
}

FullSpaceObservationResidual full_space_observation_jvp(
    const FullSpaceObservation& observation,
    const BoundaryVector& delta_reconstructed,
    const BoundaryVector& delta_state_boundary)
{
    validate(observation);

    FullSpaceObservationResidual result;

    for (std::size_t c = 0; c < boundary_channels; ++c) {
        result.constraint[c] =
            delta_reconstructed[c] -
            delta_state_boundary[c];

        if (observation.included) {
            result.measurement[c] =
                delta_reconstructed[c] /
                observation.scales[c];
        }
    }

    return result;
}

FullSpaceObservationAdjoint full_space_observation_vjp(
    const FullSpaceObservation& observation,
    const FullSpaceObservationResidual& dual)
{
    validate(observation);

    FullSpaceObservationAdjoint result;

    for (std::size_t c = 0; c < boundary_channels; ++c) {
        result.reconstructed[c] =
            dual.constraint[c];

        if (observation.included)
            result.reconstructed[c] +=
                dual.measurement[c] /
                observation.scales[c];

        result.state_boundary[c] =
            -dual.constraint[c];
    }

    return result;
}

double full_space_measurement_cost(
    const FullSpaceObservation& observation)
{
    const auto residual =
        evaluate_full_space_observation(observation);

    double cost = 0.0;

    for (double value : residual.measurement)
        cost += 0.5 * value * value;

    return cost;
}

} // namespace tpsf
