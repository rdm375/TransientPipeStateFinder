#include "tpsf/full_space_measurement_objective.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace tpsf {

MeasurementObjectiveResult
evaluate_full_space_measurement_objective(
    const std::vector<BoundaryVector>& reconstructed,
    const std::vector<BoundaryMeasurement>& measurements)
{
    if (reconstructed.size() != measurements.size()) {
        throw std::invalid_argument(
            "Measurement and reconstruction time counts differ");
    }

    MeasurementObjectiveResult result;

    result.gradient.resize(reconstructed.size());

    // Accumulate squared normalized residuals.
    double sum_squares = 0.0;
    double maximum_error = 0.0;

    for (std::size_t k = 0; k < measurements.size(); ++k) {
        for (std::size_t c = 0; c < 4; ++c) {
            const auto& measurement = measurements[k];

            if (!measurement.available[c]) {
                result.gradient[k][c] = 0.0;
                continue;
            }

            const double y = measurement.observed[c];
            const double sigma = measurement.sigma[c];
            const double b = reconstructed[k][c];

            if (!std::isfinite(y) ||
                !std::isfinite(sigma) ||
                !std::isfinite(b) ||
                !(sigma > 0.0)) {
                throw std::invalid_argument(
                    "Invalid available boundary measurement");
            }

            const double residual = b-y;
            const double normalized = residual/sigma;
            const double gradient = normalized/sigma;

            if (!std::isfinite(normalized) ||
                !std::isfinite(gradient)) {
                throw std::invalid_argument(
                    "Nonfinite normalized measurement residual "
                    "or gradient");
            }

            const double squared = normalized*normalized;

            if (!std::isfinite(squared) ||
                !std::isfinite(sum_squares+squared)) {
                throw std::invalid_argument(
                    "Measurement objective overflow");
            }

            sum_squares += squared;

            maximum_error = std::max(
                maximum_error, std::abs(normalized));

            result.gradient[k][c] = gradient;
            ++result.measurement_count;
        }
    }

    result.value = 0.5*sum_squares;

    if (result.measurement_count != 0) {
        result.normalized_rms = std::sqrt(
            sum_squares /
            static_cast<double>(result.measurement_count));

        result.maximum_normalized_error = maximum_error;
    }

    return result;
}

} // namespace tpsf
