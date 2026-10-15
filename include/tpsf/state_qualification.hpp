#pragma once

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace tpsf {

struct StateQualificationSample {
    double analysis_time_s = 0.0;

    bool reconstruction_converged = false;
    bool restart_succeeded = false;

    double pressure_error_pa = 0.0;
    double flow_error_kg_s = 0.0;

    double normalized_error = 0.0;
    bool acceptable = false;
};

struct StateQualificationOptions {
    double pressure_tolerance_pa = 1.0e4;
    double flow_tolerance_kg_s = 0.1;
    double restart_horizon_s = 4.0;
};

struct StateQualificationResult {
    StateQualificationOptions options;
    std::vector<StateQualificationSample> samples;

    std::optional<double> qualification_time_s;
    std::optional<std::size_t> qualification_index;

    bool qualified = false;
};

StateQualificationResult assess_state_qualification(
    std::vector<StateQualificationSample> samples,
    const StateQualificationOptions& options);

std::string state_qualification_json(
    const StateQualificationResult& result);

std::string state_qualification_markdown(
    const StateQualificationResult& result);

void write_state_qualification(
    const StateQualificationResult& result,
    const std::filesystem::path& directory);

} // namespace tpsf
