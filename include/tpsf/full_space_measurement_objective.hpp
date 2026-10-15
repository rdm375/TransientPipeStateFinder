#pragma once

#include "tpsf/physical_boundary_channels.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace tpsf {

struct BoundaryMeasurement {
    BoundaryVector observed{};
    BoundaryVector sigma{};
    std::array<bool, 4> available{};
};

struct MeasurementObjectiveResult {
    double value = 0.0;

    // Exact derivative with respect to reconstructed
    // boundary channels at each observation time.
    std::vector<BoundaryVector> gradient;

    std::size_t measurement_count = 0;

    // Undefined when no measurements are available.
    std::optional<double> normalized_rms;
    std::optional<double> maximum_normalized_error;
};

// Evaluates the unnormalized weighted least-squares
// objective and its exact reconstructed-channel gradient.
//
// Missing channels are ignored completely, including
// their observed and sigma fields.
//
// The two input arrays must have identical time counts.
MeasurementObjectiveResult
evaluate_full_space_measurement_objective(
    const std::vector<BoundaryVector>& reconstructed,
    const std::vector<BoundaryMeasurement>& measurements);

} // namespace tpsf
