#include "tpsf/sensitivity_decay.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tpsf {
namespace {

AttenuationResult persistent_attenuation(
    const std::vector<double>& time,
    const std::vector<double>& sensitivity,
    double threshold)
{
    AttenuationResult result;
    result.threshold = threshold;

    // Scan backward to find the beginning of the final
    // suffix that remains below the threshold.
    //
    // A temporary crossing is not sufficient.
    std::size_t suffix_start = sensitivity.size();

    for (std::size_t k = sensitivity.size(); k-- > 0;) {
        if (sensitivity[k] > threshold)
            break;

        suffix_start = k;
    }

    if (suffix_start < sensitivity.size()) {
        result.status = AttenuationStatus::attained;
        result.index = suffix_start;
        result.time = time[suffix_start];
    }

    return result;
}

} // namespace

SensitivityDecayHistory assess_sensitivity_decay(
    const std::vector<double>& time,
    const std::vector<double>& sensitivity,
    const SensitivityDecayOptions& options)
{
    if (time.empty() || time.size() != sensitivity.size())
        throw std::invalid_argument(
            "Invalid sensitivity history dimensions");

    if (!std::isfinite(options.absolute_threshold) ||
        options.absolute_threshold < 0.0 ||
        !std::isfinite(options.relative_threshold) ||
        options.relative_threshold < 0.0)
        throw std::invalid_argument(
            "Invalid sensitivity decay thresholds");

    for (std::size_t k = 0; k < time.size(); ++k) {
        if (!std::isfinite(time[k]) ||
            !std::isfinite(sensitivity[k]) ||
            sensitivity[k] < 0.0)
            throw std::invalid_argument(
                "Nonfinite or negative sensitivity history");

        if (k > 0 && !(time[k] > time[k-1]))
            throw std::invalid_argument(
                "Forecast times must be strictly increasing");
    }

    SensitivityDecayHistory result;

    result.initial = sensitivity.front();
    result.terminal = sensitivity.back();

    const auto peak = std::max_element(
        sensitivity.begin(), sensitivity.end());

    result.maximum = *peak;

    result.peak_index = static_cast<std::size_t>(
        peak - sensitivity.begin());

    result.peak_time = time[result.peak_index];

    result.delayed_peak = result.peak_index > 0;

    for (std::size_t k = 1; k < sensitivity.size(); ++k) {
        if (sensitivity[k] > sensitivity[k-1]) {
            result.nonmonotonic = true;
            break;
        }
    }

    result.absolute = persistent_attenuation(
        time, sensitivity, options.absolute_threshold);

    if (result.initial > 0.0) {
        result.relative = persistent_attenuation(
            time,
            sensitivity,
            options.relative_threshold * result.initial);
    } else {
        result.relative.status =
            AttenuationStatus::undefined_reference;
    }

    // Peak-relative attenuation uses the maximum
    // sensitivity over the simulated forecast horizon.
    //
    // This remains meaningful for delayed responses
    // whose initial sensitivity is zero.
    if (result.maximum > 0.0) {
        result.peak_relative = persistent_attenuation(
            time,
            sensitivity,
            options.relative_threshold * result.maximum);
    } else {
        result.peak_relative.status =
            AttenuationStatus::undefined_reference;
    }

    return result;
}

ForecastSensitivityDecay assess_forecast_sensitivity_decay(
    const ForecastDirectionalUncertainty& forecast,
    const SensitivityDecayOptions& options)
{
    ForecastSensitivityDecay result;

    result.measurement = assess_sensitivity_decay(
        forecast.time,
        forecast.measurement_unit_radius,
        options);

    result.unresolved = assess_sensitivity_decay(
        forecast.time,
        forecast.unresolved_unit_radius,
        options);

    return result;
}

} // namespace tpsf
