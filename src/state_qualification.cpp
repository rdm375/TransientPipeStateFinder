#include "tpsf/state_qualification.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace tpsf {

StateQualificationResult assess_state_qualification(
    std::vector<StateQualificationSample> samples,
    const StateQualificationOptions& options)
{
    if (!(std::isfinite(options.pressure_tolerance_pa) &&
          options.pressure_tolerance_pa > 0.0 &&
          std::isfinite(options.flow_tolerance_kg_s) &&
          options.flow_tolerance_kg_s > 0.0 &&
          std::isfinite(options.restart_horizon_s) &&
          options.restart_horizon_s > 0.0))
        throw std::invalid_argument(
            "Invalid state qualification tolerances or horizon");

    if (samples.empty())
        throw std::invalid_argument(
            "State qualification requires candidate states");

    StateQualificationResult result;
    result.options = options;
    result.samples = std::move(samples);

    double previous_time =
        -std::numeric_limits<double>::infinity();

    for (auto& sample : result.samples) {
        if (!std::isfinite(sample.analysis_time_s) ||
            sample.analysis_time_s <= previous_time)
            throw std::invalid_argument(
                "Qualification times must be finite and increasing");

        previous_time = sample.analysis_time_s;

        sample.acceptable = false;
        sample.normalized_error =
            std::numeric_limits<double>::infinity();

        if (!sample.reconstruction_converged ||
            !sample.restart_succeeded)
            continue;

        if (!std::isfinite(sample.pressure_error_pa) ||
            !std::isfinite(sample.flow_error_kg_s) ||
            sample.pressure_error_pa < 0.0 ||
            sample.flow_error_kg_s < 0.0)
            continue;

        sample.normalized_error = std::max(
            sample.pressure_error_pa /
                options.pressure_tolerance_pa,
            sample.flow_error_kg_s /
                options.flow_tolerance_kg_s);

        sample.acceptable =
            sample.normalized_error <= 1.0;
    }

    // A candidate qualifies persistently only when
    // it and every later tested candidate are acceptable.
    bool acceptable_suffix = true;

    for (std::size_t k = result.samples.size();
         k-- > 0;) {

        acceptable_suffix =
            acceptable_suffix &&
            result.samples[k].acceptable;

        if (acceptable_suffix) {
            result.qualification_index = k;
            result.qualification_time_s =
                result.samples[k].analysis_time_s;
            result.qualified = true;
        }
    }

    return result;
}

std::string state_qualification_json(
    const StateQualificationResult& result)
{
    std::ostringstream out;
    out << std::setprecision(17);

    out << "{\n";
    out << "  \"schema\": \"tpsf.state-qualification.v1\",\n";
    out << "  \"primary_engineering_result\": true,\n";
    out << "  \"method\": \"causal reconstruction and nonlinear restart\",\n";
    out << "  \"reference_basis\": \"synthetic truth\",\n";
    out << "  \"qualification_status\": \""
        << (result.qualified ? "qualified" : "not_qualified")
        << "\",\n";

    out << "  \"qualification_time_s\": ";

    if (result.qualification_time_s)
        out << *result.qualification_time_s;
    else
        out << "null";

    out << ",\n";

    out << "  \"pressure_tolerance_pa\": "
        << result.options.pressure_tolerance_pa << ",\n";
    out << "  \"flow_tolerance_kg_s\": "
        << result.options.flow_tolerance_kg_s << ",\n";
    out << "  \"restart_horizon_s\": "
        << result.options.restart_horizon_s << ",\n";

    out << "  \"samples\": [\n";

    for (std::size_t k = 0;
         k < result.samples.size(); ++k) {

        const auto& s = result.samples[k];

        out << "    {\n";
        out << "      \"analysis_time_s\": "
            << s.analysis_time_s << ",\n";
        out << "      \"reconstruction_converged\": "
            << (s.reconstruction_converged ? "true" : "false")
            << ",\n";
        out << "      \"restart_succeeded\": "
            << (s.restart_succeeded ? "true" : "false")
            << ",\n";

        out << "      \"pressure_error_pa\": ";

        if (s.restart_succeeded &&
            std::isfinite(s.pressure_error_pa))
            out << s.pressure_error_pa;
        else
            out << "null";

        out << ",\n";
        out << "      \"flow_error_kg_s\": ";

        if (s.restart_succeeded &&
            std::isfinite(s.flow_error_kg_s))
            out << s.flow_error_kg_s;
        else
            out << "null";

        out << ",\n";
        out << "      \"normalized_error\": ";

        if (std::isfinite(s.normalized_error))
            out << s.normalized_error;
        else
            out << "null";

        out << ",\n";
        out << "      \"acceptable\": "
            << (s.acceptable ? "true" : "false")
            << "\n";
        out << "    }";

        if (k + 1 < result.samples.size())
            out << ",";

        out << "\n";
    }

    out << "  ],\n";
    out << "  \"limitations\": [\n";
    out << "    \"Synthetic reference, not field validation\",\n";
    out << "    \"Qualification applies to tested analysis times\",\n";
    out << "    \"Qualification applies to the specified restart horizon\",\n";
    out << "    \"Qualification applies to the tested future controls\",\n";
    out << "    \"No guarantee for arbitrary future operation\"\n";
    out << "  ]\n";
    out << "}\n";

    return out.str();
}

std::string state_qualification_markdown(
    const StateQualificationResult& result)
{
    std::ostringstream out;
    out << std::setprecision(12);

    out << "# State Qualification Assessment\n\n";
    out << "## PRIMARY ENGINEERING RESULT\n\n";

    out << "**Qualification status:** "
        << (result.qualified
            ? "QUALIFIED"
            : "NOT QUALIFIED")
        << "\n\n";

    out << "**Earliest persistently qualified analysis time:** ";

    if (result.qualification_time_s)
        out << *result.qualification_time_s << " s\n\n";
    else
        out << "Not attained\n\n";

    out << "A qualified state can initialize a new nonlinear "
           "simulation satisfying the specified pressure "
           "and mass-flow error tolerances over the tested "
           "restart horizon, relative to synthetic truth.\n\n";

    out << "### Qualification requirements\n\n";
    out << "- Pressure tolerance: "
        << result.options.pressure_tolerance_pa
        << " Pa\n";
    out << "- Mass-flow tolerance: "
        << result.options.flow_tolerance_kg_s
        << " kg/s\n";
    out << "- Nonlinear restart horizon: "
        << result.options.restart_horizon_s
        << " s\n";
    out << "- Reference: synthetic truth\n";
    out << "- Reconstruction: causal measurement prefixes\n";
    out << "- Acceptance: persistent across later tested "
           "analysis times\n\n";

    out << "### Qualification history\n\n";

    out << "| Analysis time (s) | Pressure error (Pa) | "
           "Flow error (kg/s) | Normalized error | "
           "Converged | Acceptable |\n";
    out << "|---:|---:|---:|---:|---|---|\n";

    for (const auto& s : result.samples) {
        out << "| " << s.analysis_time_s << " | ";

        if (s.restart_succeeded) {
            out << s.pressure_error_pa
                << " | " << s.flow_error_kg_s;
        } else {
            out << "unavailable | unavailable";
        }

        out << " | ";

        if (std::isfinite(s.normalized_error))
            out << s.normalized_error;
        else
            out << "unavailable";

        out << " | "
            << (s.reconstruction_converged ? "yes" : "no")
            << " | "
            << (s.acceptable ? "yes" : "no")
            << " |\n";
    }

    out << "\n### Interpretation\n\n";
    out << "The qualification time is the first sampled "
           "analysis time for which every remaining tested "
           "candidate passes the nonlinear restart test.\n\n";
    out << "The result is conditional on the stated tolerances, "
           "future boundary conditions, restart horizon, "
           "and synthetic reference. It is not an "
           "unrestricted guarantee of future accuracy.\n";

    return out.str();
}

void write_state_qualification(
    const StateQualificationResult& result,
    const std::filesystem::path& directory)
{
    const auto json = state_qualification_json(result);
    const auto markdown =
        state_qualification_markdown(result);

    std::filesystem::create_directories(directory);

    std::ofstream json_file(
        directory / "state_qualification.json");

    std::ofstream markdown_file(
        directory / "state_qualification.md");

    if (!json_file || !markdown_file)
        throw std::runtime_error(
            "Cannot create qualification reports");

    json_file << json;
    markdown_file << markdown;

    if (!json_file || !markdown_file)
        throw std::runtime_error(
            "Cannot write qualification reports");
}

} // namespace tpsf
