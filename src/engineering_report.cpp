#include "tpsf/engineering_report.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace tpsf {
namespace {

std::string escape_json(const std::string& input)
{
    std::string result;

    for (unsigned char c : input) {
        switch (c) {
        case '"': result += "\\\""; break;
        case '\\': result += "\\\\"; break;
        case '\n': result += "\\n"; break;
        case '\r': result += "\\r"; break;
        case '\t': result += "\\t"; break;
        default:
            if (c < 0x20) {
                constexpr char hex[] = "0123456789abcdef";
                result += "\\u00";
                result += hex[c >> 4];
                result += hex[c & 15];
            } else {
                result += static_cast<char>(c);
            }
        }
    }

    return result;
}

std::string quoted(const std::string& value)
{
    return "\"" + escape_json(value) + "\"";
}

std::string number(double value)
{
    if (!std::isfinite(value))
        return "null";

    std::ostringstream stream;
    stream << std::setprecision(17) << value;
    return stream.str();
}

std::string boolean(bool value)
{
    return value ? "true" : "false";
}

std::string optional_number(
    const std::optional<double>& value)
{
    return value ? number(*value) : "null";
}

std::string classification(
    DecisionClassification value)
{
    switch (value) {
    case DecisionClassification::certified_satisfactory:
        return "satisfactory";
    case DecisionClassification::certified_violation:
        return "violation";
    case DecisionClassification::indeterminate:
        return "indeterminate";
    }

    return "unknown";
}

std::string markdown_classification(
    DecisionClassification value)
{
    return classification(value);
}

std::string optional_time(
    const std::optional<double>& value)
{
    return value ? number(*value) + " s" : "not attained";
}

std::string markdown_number(double value)
{
    return std::isfinite(value)
        ? number(value)
        : (std::isnan(value)
            ? "undefined"
            : (value < 0 ? "-infinity" : "+infinity"));
}

void validate(const EngineeringReport& report)
{
    const auto& s = report.scenario;
    const auto& d = report.decision.directional;
    const auto& points = report.decision.decision.points;

    if (s.cells <= 0 ||
        s.reconstruction_steps < 0 ||
        s.forecast_steps < 1 ||
        !std::isfinite(s.dt) ||
        s.dt <= 0.0 ||
        !std::isfinite(s.theta) ||
        s.theta < 0.0 ||
        s.theta > 1.0 ||
        !std::isfinite(s.minimum_outlet_pressure_pa) ||
        !std::isfinite(s.measurement_radius) ||
        s.measurement_radius < 0.0 ||
        !std::isfinite(s.unresolved_radius) ||
        s.unresolved_radius < 0.0)
        throw std::invalid_argument(
            "Invalid engineering scenario");

    if (points.empty() ||
        points.size() != d.time.size() ||
        points.size() != d.nominal_outlet_pressure.size())
        throw std::invalid_argument(
            "Decision/forecast history size mismatch");

    if (report.quality.information_singular_values.size() < 0)
        throw std::invalid_argument("Invalid information spectrum");

    for (const auto& m : s.measurements) {
        for (std::size_t j = 0; j < 4; ++j) {
            if (m.available[j] &&
                (!std::isfinite(m.sigma[j]) ||
                 m.sigma[j] <= 0.0))
                throw std::invalid_argument(
                    "Invalid available measurement sigma");
        }
    }
}

} // namespace

std::string engineering_report_json(
    const EngineeringReport& report)
{
    validate(report);

    const auto& p = report.provenance;
    const auto& s = report.scenario;
    const auto& q = report.quality;
    const auto& d = report.decision;

    const auto& singular = q.information_singular_values;

    double largest = 0.0;
    double smallest = std::numeric_limits<double>::infinity();

    for (Eigen::Index i = 0; i < singular.size(); ++i) {
        const double value = singular[i];

        if (!std::isfinite(value) || value < 0.0)
            throw std::invalid_argument(
                "Invalid information singular value");

        largest = std::max(largest, value);

        if (value > 0.0 &&
            value > q.relative_rank_tolerance * singular[0])
            smallest = std::min(smallest, value);
    }

    const bool finite_condition =
        q.feasible_dimension > 0 &&
        q.measurement_rank == q.feasible_dimension &&
        largest > 0.0 &&
        std::isfinite(smallest);

    const double condition =
        finite_condition ? largest / smallest
                         : std::numeric_limits<double>::infinity();

    std::ostringstream out;
    out << std::setprecision(17);

    out << "{\n";
    out << "  \"schema\": \"tpsf.engineering-assessment.v1\",\n";

    out << "  \"provenance\": {\n";
    out << "    \"statefinder_revision\": "
        << quoted(p.statefinder_revision) << ",\n";
    out << "    \"solver_revision\": "
        << quoted(p.solver_revision) << ",\n";
    out << "    \"validated_solver_baseline\": "
        << quoted(p.validated_solver_baseline) << ",\n";
    out << "    \"statefinder_dirty\": "
        << boolean(p.statefinder_dirty) << ",\n";
    out << "    \"solver_dirty\": "
        << boolean(p.solver_dirty) << ",\n";
    out << "    \"build_type\": "
        << quoted(p.build_type) << ",\n";
    out << "    \"compiler\": "
        << quoted(p.compiler) << "\n";
    out << "  },\n";

    out << "  \"scenario\": {\n";
    out << "    \"cells\": " << s.cells << ",\n";
    out << "    \"reconstruction_steps\": "
        << s.reconstruction_steps << ",\n";
    out << "    \"forecast_steps\": "
        << s.forecast_steps << ",\n";
    out << "    \"dt_s\": " << number(s.dt) << ",\n";
    out << "    \"theta\": " << number(s.theta) << ",\n";
    out << "    \"minimum_outlet_pressure_pa\": "
        << number(s.minimum_outlet_pressure_pa) << ",\n";
    out << "    \"measurement_radius\": "
        << number(s.measurement_radius) << ",\n";
    out << "    \"bound_unresolved\": "
        << boolean(s.bound_unresolved) << ",\n";
    out << "    \"unresolved_radius\": "
        << number(s.unresolved_radius) << ",\n";

    out << "    \"measurements\": [\n";

    for (std::size_t k = 0; k < s.measurements.size(); ++k) {
        const auto& m = s.measurements[k];

        out << "      {\"available\": [";

        for (std::size_t j = 0; j < 4; ++j) {
            if (j) out << ", ";
            out << boolean(m.available[j]);
        }

        out << "], \"sigma\": [";

        for (std::size_t j = 0; j < 4; ++j) {
            if (j) out << ", ";
            out << number(m.sigma[j]);
        }

        out << "]}";

        if (k + 1 != s.measurements.size())
            out << ",";

        out << "\n";
    }

    out << "    ]\n";
    out << "  },\n";

    out << "  \"reconstruction\": {\n";
    out << "    \"converged\": "
        << boolean(q.converged) << ",\n";
    out << "    \"objective\": "
        << number(q.objective) << ",\n";
    out << "    \"constraint_l2\": "
        << number(q.constraint_l2) << ",\n";
    out << "    \"measurement_l2\": "
        << number(q.measurement_l2) << ",\n";
    out << "    \"projected_gradient\": "
        << number(q.projected_gradient) << ",\n";
    out << "    \"feasible_dimension\": "
        << q.feasible_dimension << ",\n";
    out << "    \"measurement_rank\": "
        << q.measurement_rank << ",\n";
    out << "    \"unresolved_directions\": "
        << q.unresolved_directions << ",\n";
    out << "    \"locally_identifiable\": "
        << boolean(q.locally_identifiable) << ",\n";
    out << "    \"rank_tolerance\": "
        << number(q.relative_rank_tolerance) << ",\n";
    out << "    \"information_condition_number\": "
        << number(condition) << ",\n";
    out << "    \"information_condition_finite\": "
        << boolean(finite_condition) << ",\n";
    out << "    \"information_singular_values\": [";

    for (Eigen::Index i = 0; i < singular.size(); ++i) {
        if (i) out << ", ";
        out << number(singular[i]);
    }

    out << "]\n";
    out << "  },\n";

    out << "  \"uncertainty\": {\n";
    out << "    \"model\": \"local-linearized deterministic\",\n";
    out << "    \"statistical_confidence_claim\": false,\n";
    out << "    \"measurement_coordinate_radius\": "
        << number(d.measurement_coordinate_radius) << ",\n";
    out << "    \"unresolved_bounded\": "
        << boolean(d.unresolved_bounded) << "\n";
    out << "  },\n";

    out << "  \"decision\": {\n";
    out << "    \"classification\": "
        << quoted(classification(d.decision.classification))
        << ",\n";
    out << "    \"has_reliable_suffix\": "
        << boolean(d.decision.has_reliable_suffix) << ",\n";
    out << "    \"reliable_suffix_start\": ";

    if (d.decision.has_reliable_suffix)
        out << d.decision.reliable_suffix_start;
    else
        out << "null";

    out << ",\n";
    out << "    \"points\": [\n";

    for (std::size_t k = 0;
         k < d.decision.points.size(); ++k) {

        const auto& point = d.decision.points[k];

        out << "      {\n";
        out << "        \"time_s\": "
            << number(d.directional.time[k]) << ",\n";
        out << "        \"classification\": "
            << quoted(classification(point.classification))
            << ",\n";
        out << "        \"nominal_pa\": "
            << number(point.nominal) << ",\n";
        out << "        \"limit_pa\": "
            << number(point.limit) << ",\n";
        out << "        \"lower_bound_pa\": "
            << number(point.lower_bound) << ",\n";
        out << "        \"upper_bound_pa\": "
            << number(point.upper_bound) << ",\n";
        out << "        \"uncertainty_radius_pa\": "
            << number(point.uncertainty_radius) << ",\n";
        out << "        \"worst_case_margin_pa\": "
            << number(point.worst_case_margin) << ",\n";
        out << "        \"finite_bounds\": "
            << boolean(point.finite_bounds) << "\n";
        out << "      }";

        if (k + 1 != d.decision.points.size())
            out << ",";

        out << "\n";
    }

    out << "    ]\n";
    out << "  },\n";

    out << "  \"memory\": {\n";
    out << "    \"relative_attenuation_fraction\": "
        << number(report.attenuation_fraction) << ",\n";
    out << "    \"state_persistent_time_s\": "
        << optional_number(report.state_attenuation_time)
        << ",\n";
    out << "    \"boundary_persistent_time_s\": "
        << optional_number(report.boundary_attenuation_time)
        << ",\n";
    out << "    \"outlet_persistent_time_s\": "
        << optional_number(report.outlet_attenuation_time)
        << ",\n";
    out << "    \"maximum_tangent_residual\": "
        << number(report.maximum_memory_tangent_residual)
        << ",\n";
    out << "    \"maximum_mass_defect\": "
        << number(report.maximum_memory_mass_defect)
        << "\n";
    out << "  },\n";

    out << "  \"limitations\": [\n";
    out << "    \"Synthetic demonstration\",\n";
    out << "    \"Local-linearized uncertainty only\",\n";
    out << "    \"No nonlinear decision certificate\",\n";
    out << "    \"Finite sampled forecast horizon\",\n";
    out << "    \"No statistical confidence interpretation\"\n";
    out << "  ]\n";
    out << "}\n";

    return out.str();
}

std::string engineering_report_markdown(
    const EngineeringReport& report)
{
    validate(report);

    const auto& p = report.provenance;
    const auto& s = report.scenario;
    const auto& q = report.quality;
    const auto& d = report.decision;

    std::ostringstream out;
    out << std::setprecision(12);

    out << "# TransientPipeStateFinder Engineering Assessment\n\n";
    out << "**Status:** Research demonstration; "
           "not an operational certificate.\n\n";

    out << "## Provenance\n\n";
    out << "- Statefinder revision: `"
        << p.statefinder_revision << "`\n";
    out << "- Solver revision: `"
        << p.solver_revision << "`\n";
    out << "- Validated solver baseline: `"
        << p.validated_solver_baseline << "`\n";
    out << "- Statefinder working tree modified: "
        << (p.statefinder_dirty ? "yes/unknown" : "no")
        << "\n";
    out << "- Solver working tree modified: "
        << (p.solver_dirty ? "yes/unknown" : "no")
        << "\n";
    out << "- Build type: " << p.build_type << "\n";
    out << "- Compiler: " << p.compiler << "\n\n";

    out << "## Scenario\n\n";
    out << "- Cells: " << s.cells << "\n";
    out << "- Reconstruction steps: "
        << s.reconstruction_steps << "\n";
    out << "- Forecast steps: "
        << s.forecast_steps << "\n";
    out << "- Time step: " << s.dt << " s\n";
    out << "- Theta-method weighting: "
        << s.theta << "\n";
    out << "- Minimum outlet pressure: "
        << s.minimum_outlet_pressure_pa << " Pa\n";
    out << "- Standardized measurement radius: "
        << s.measurement_radius << "\n";
    out << "- Unresolved directions bounded: "
        << (s.bound_unresolved ? "yes" : "no") << "\n\n";

    out << "### Measurement channels\n\n";
    out << "Channel order: inlet pressure, inlet flow, "
           "outlet pressure, outlet flow.\n\n";
    out << "| Level | Available | Sigma |\n";
    out << "|---:|---|---|\n";

    for (std::size_t k = 0; k < s.measurements.size(); ++k) {
        const auto& m = s.measurements[k];

        out << "| " << k << " | ";

        for (std::size_t j = 0; j < 4; ++j)
            out << (m.available[j] ? "1" : "0");

        out << " | ";

        for (std::size_t j = 0; j < 4; ++j) {
            if (j) out << ", ";
            out << m.sigma[j];
        }

        out << " |\n";
    }

    out << "\n## Reconstruction and identifiability\n\n";
    out << "- Converged: "
        << (q.converged ? "yes" : "no") << "\n";
    out << "- Constraint L2: "
        << q.constraint_l2 << "\n";
    out << "- Measurement L2: "
        << q.measurement_l2 << "\n";
    out << "- Feasible dimension: "
        << q.feasible_dimension << "\n";
    out << "- Measurement rank: "
        << q.measurement_rank << "\n";
    out << "- Unresolved directions: "
        << q.unresolved_directions << "\n";
    out << "- Locally identifiable: "
        << (q.locally_identifiable ? "yes" : "no")
        << "\n";
    out << "- Relative rank tolerance: "
        << q.relative_rank_tolerance << "\n";

    out << "\nInformation singular values:\n\n";

    for (Eigen::Index i = 0;
         i < q.information_singular_values.size(); ++i) {

        out << "- " << q.information_singular_values[i]
            << "\n";
    }

    out << "\n## Decision reliability\n\n";
    out << "- Overall classification: **"
        << markdown_classification(
            d.decision.classification)
        << "**\n";
    out << "- Reliable suffix: "
        << (d.decision.has_reliable_suffix ? "yes" : "no")
        << "\n\n";

    out << "| Time (s) | Nominal (Pa) | Lower (Pa) | "
           "Upper (Pa) | Worst margin (Pa) | Classification |\n";
    out << "|---:|---:|---:|---:|---:|---|\n";

    for (std::size_t k = 0;
         k < d.decision.points.size(); ++k) {

        const auto& point = d.decision.points[k];

        out << "| " << d.directional.time[k]
            << " | " << markdown_number(point.nominal)
            << " | " << markdown_number(point.lower_bound)
            << " | " << markdown_number(point.upper_bound)
            << " | " << markdown_number(point.worst_case_margin)
            << " | " << classification(point.classification)
            << " |\n";
    }

    out << "\n## Initial-state memory\n\n";
    out << "- Relative attenuation threshold: "
        << report.attenuation_fraction << "\n";
    out << "- Full-state persistent attenuation: "
        << optional_time(report.state_attenuation_time)
        << "\n";
    out << "- Boundary persistent attenuation: "
        << optional_time(report.boundary_attenuation_time)
        << "\n";
    out << "- Outlet-pressure persistent attenuation: "
        << optional_time(report.outlet_attenuation_time)
        << "\n";
    out << "- Maximum tangent equation residual: "
        << report.maximum_memory_tangent_residual
        << "\n";
    out << "- Maximum tangent mass defect: "
        << report.maximum_memory_mass_defect
        << "\n\n";

    out << "## Interpretation and limitations\n\n";
    out << "Numerical convergence, local identifiability, "
           "and engineering decision reliability are distinct.\n\n";
    out << "Uncertainty is represented by a deterministic "
           "local-linearized model. The measurement radius "
           "is not a statistical confidence level.\n\n";
    out << "Unresolved directions must not be assumed bounded "
           "unless explicitly specified.\n\n";
    out << "Attenuation times apply only to the sampled "
           "forecast horizon. They do not establish "
           "subsequent behavior.\n\n";
    out << "This assessment is not a nonlinear certificate "
           "or an operational authorization.\n";

    return out.str();
}

void write_engineering_report(
    const EngineeringReport& report,
    const std::filesystem::path& directory)
{
    // Render both before touching output files.
    const auto json = engineering_report_json(report);
    const auto markdown = engineering_report_markdown(report);

    std::filesystem::create_directories(directory);

    const auto json_path =
        directory / "engineering_assessment.json";

    const auto markdown_path =
        directory / "engineering_assessment.md";

    std::ofstream json_file(json_path);
    std::ofstream markdown_file(markdown_path);

    if (!json_file || !markdown_file)
        throw std::runtime_error(
            "Unable to create engineering assessment files");

    json_file << json;
    markdown_file << markdown;

    if (!json_file || !markdown_file)
        throw std::runtime_error(
            "Unable to write engineering assessment files");
}

} // namespace tpsf
