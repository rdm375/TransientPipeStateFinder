#pragma once

#include "tpsf/forecast_directional_uncertainty.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace tpsf {

enum class AttenuationStatus {
    attained,
    not_attained,
    undefined_reference
};

struct AttenuationResult {
    AttenuationStatus status = AttenuationStatus::not_attained;
    std::optional<std::size_t> index;
    std::optional<double> time;
    double threshold = 0.0;
};

struct SensitivityDecayOptions {
    // Absolute thresholds have the units of the quantity
    // being analyzed (Pa for outlet-pressure uncertainty).
    double absolute_threshold = 0.0;

    // Relative threshold is a fraction of the initial value.
    double relative_threshold = 0.05;
};

struct SensitivityDecayHistory {
    double initial = 0.0;
    double maximum = 0.0;
    double terminal = 0.0;

    std::size_t peak_index = 0;
    double peak_time = 0.0;

    bool nonmonotonic = false;
    bool delayed_peak = false;

    AttenuationResult absolute;

    // Relative to the sensitivity at the initial time.
    AttenuationResult relative;

    // Relative to the maximum sensitivity over the
    // assessed forecast horizon.
    AttenuationResult peak_relative;
};

struct ForecastSensitivityDecay {
    SensitivityDecayHistory measurement;
    SensitivityDecayHistory unresolved;
};

SensitivityDecayHistory assess_sensitivity_decay(
    const std::vector<double>& time,
    const std::vector<double>& sensitivity,
    const SensitivityDecayOptions& options = {});

ForecastSensitivityDecay assess_forecast_sensitivity_decay(
    const ForecastDirectionalUncertainty& forecast,
    const SensitivityDecayOptions& options = {});

} // namespace tpsf
