#pragma once

#include "tpsf/state_qualification.hpp"
#include "tpsf/statefinding_quality.hpp"
#include "tpsf/pipeline_decision.hpp"

#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace tpsf {

struct EngineeringProvenance {
    std::string statefinder_revision = "unknown";
    std::string solver_revision = "unknown";

    std::string validated_solver_baseline =
        "670ac493e4461e678ca66f2254dc3d402df3f312";

    bool statefinder_dirty = true;
    bool solver_dirty = true;

    std::string build_type = "unknown";
    std::string compiler = "unknown";
};

struct EngineeringMeasurement {
    std::array<bool, 4> available{};
    std::array<double, 4> sigma{};
};

struct EngineeringScenario {
    int cells = 0;
    int reconstruction_steps = 0;
    int forecast_steps = 0;

    double dt = 0.0;
    double theta = 0.0;

    double minimum_outlet_pressure_pa = 0.0;
    double measurement_radius = 1.0;

    bool bound_unresolved = false;
    double unresolved_radius = 0.0;

    std::vector<EngineeringMeasurement> measurements;
};

struct EngineeringReport {
    std::optional<StateQualificationResult> qualification;
    EngineeringProvenance provenance;
    EngineeringScenario scenario;

    StatefindingQuality quality;
    PipelineDecisionResult decision;

    double maximum_memory_tangent_residual = 0.0;
    double maximum_memory_mass_defect = 0.0;

    std::optional<double> state_attenuation_time;
    std::optional<double> boundary_attenuation_time;
    std::optional<double> outlet_attenuation_time;

    double attenuation_fraction = 0.05;
};

std::string engineering_report_json(
    const EngineeringReport& report);

std::string engineering_report_markdown(
    const EngineeringReport& report);

void write_engineering_report(
    const EngineeringReport& report,
    const std::filesystem::path& directory);

} // namespace tpsf
