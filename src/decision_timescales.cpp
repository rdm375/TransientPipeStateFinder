#include "tpsf/decision_timescales.hpp"

#include <cmath>
#include <stdexcept>

namespace tpsf {
namespace {

void validate_time(const std::vector<double>& time)
{
    if (time.empty())
        throw std::invalid_argument("Empty forecast time grid");

    for (std::size_t k = 0; k < time.size(); ++k) {
        if (!std::isfinite(time[k]))
            throw std::invalid_argument(
                "Nonfinite forecast time");

        if (k > 0 && !(time[k] > time[k - 1]))
            throw std::invalid_argument(
                "Forecast times must be strictly increasing");
    }
}

PersistentDecisionTime persistent_suffix(
    const std::vector<double>& time,
    const ForecastDecisionResult& forecast,
    DecisionClassification target)
{
    PersistentDecisionTime result;

    for (std::size_t k = forecast.points.size(); k-- > 0;) {
        if (forecast.points[k].classification != target)
            break;

        result.index = k;
        result.time = time[k];
    }

    return result;
}

std::optional<double> lag(
    const PersistentDecisionTime& resolution,
    const AttenuationResult& attenuation)
{
    if (!resolution.time || !attenuation.time)
        return std::nullopt;

    return *resolution.time - *attenuation.time;
}

} // namespace

DecisionTimescales assess_decision_timescales(
    const std::vector<double>& time,
    const ForecastDecisionResult& forecast)
{
    validate_time(time);

    if (time.size() != forecast.points.size())
        throw std::invalid_argument(
            "Decision/time dimension mismatch");

    DecisionTimescales result;

    result.satisfactory = persistent_suffix(
        time,
        forecast,
        DecisionClassification::certified_satisfactory);

    result.violation = persistent_suffix(
        time,
        forecast,
        DecisionClassification::certified_violation);

    // A resolved suffix permits either satisfactory
    // or violated classifications, but no indeterminate
    // classifications.
    for (std::size_t k = forecast.points.size(); k-- > 0;) {
        if (forecast.points[k].classification ==
            DecisionClassification::indeterminate)
            break;

        result.resolved.index = k;
        result.resolved.time = time[k];
    }

    for (std::size_t k = 0; k < forecast.points.size(); ++k) {
        switch (forecast.points[k].classification) {
        case DecisionClassification::certified_satisfactory:
            ++result.satisfactory_count;
            break;

        case DecisionClassification::certified_violation:
            ++result.violation_count;
            break;

        case DecisionClassification::indeterminate:
            ++result.indeterminate_count;
            result.indeterminate_indices.push_back(k);

            if (!result.first_indeterminate_time)
                result.first_indeterminate_time = time[k];

            result.last_indeterminate_time = time[k];
            break;

        default:
            throw std::invalid_argument(
                "Unknown decision classification");
        }
    }

    return result;
}

OperationalTimescales assess_operational_timescales(
    const PipelineDecisionResult& pipeline,
    const SensitivityDecayOptions& options)
{
    OperationalTimescales result;

    result.decision = assess_decision_timescales(
        pipeline.directional.time,
        pipeline.decision);

    result.sensitivity = assess_forecast_sensitivity_decay(
        pipeline.directional,
        options);

    result.measurement_resolution_lag = lag(
        result.decision.resolved,
        result.sensitivity.measurement.peak_relative);

    result.unresolved_resolution_lag = lag(
        result.decision.resolved,
        result.sensitivity.unresolved.peak_relative);

    return result;
}

} // namespace tpsf
