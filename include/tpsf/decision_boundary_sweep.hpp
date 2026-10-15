#pragma once

#include "tpsf/pipeline_decision.hpp"

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace tpsf {

struct BoundarySweepPoint {
    double minimum_pressure_pa = 0.0;

    DecisionClassification classification =
        DecisionClassification::indeterminate;

    std::size_t satisfactory_count = 0;
    std::size_t violation_count = 0;
    std::size_t indeterminate_count = 0;

    std::optional<double> persistent_satisfactory_time;
    std::optional<double> persistent_violation_time;
};

struct DecisionBoundarySweep {
    std::vector<BoundarySweepPoint> points;

    // True when the sweep reuses an existing forecast
    // rather than rerunning nonlinear integration.
    bool reused_forecast = true;
};

DecisionBoundarySweep assess_decision_boundary_sweep(
    const PipelineDecisionResult& baseline,
    const std::vector<double>& minimum_pressures);

std::string decision_boundary_sweep_json(
    const DecisionBoundarySweep& sweep);

std::string decision_boundary_sweep_markdown(
    const DecisionBoundarySweep& sweep);

void write_decision_boundary_sweep(
    const DecisionBoundarySweep& sweep,
    const std::filesystem::path& directory);

} // namespace tpsf
