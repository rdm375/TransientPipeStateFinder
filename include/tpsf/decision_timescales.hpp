#pragma once

#include "tpsf/decision_reliability.hpp"
#include "tpsf/pipeline_decision.hpp"
#include "tpsf/sensitivity_decay.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace tpsf {

struct PersistentDecisionTime {
    std::optional<std::size_t> index;
    std::optional<double> time;
};

struct DecisionTimescales {
    PersistentDecisionTime satisfactory;
    PersistentDecisionTime violation;

    // First sample after which no indeterminate
    // classifications remain.
    PersistentDecisionTime resolved;

    std::vector<std::size_t> indeterminate_indices;

    std::optional<double> first_indeterminate_time;
    std::optional<double> last_indeterminate_time;

    std::size_t satisfactory_count = 0;
    std::size_t violation_count = 0;
    std::size_t indeterminate_count = 0;
};

struct OperationalTimescales {
    DecisionTimescales decision;
    ForecastSensitivityDecay sensitivity;

    // Difference between persistent decision resolution
    // and peak-relative attenuation, when both exist.
    std::optional<double> measurement_resolution_lag;
    std::optional<double> unresolved_resolution_lag;
};

DecisionTimescales assess_decision_timescales(
    const std::vector<double>& time,
    const ForecastDecisionResult& forecast);

OperationalTimescales assess_operational_timescales(
    const PipelineDecisionResult& pipeline,
    const SensitivityDecayOptions& options = {});

} // namespace tpsf
