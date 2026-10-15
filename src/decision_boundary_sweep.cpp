#include "tpsf/decision_boundary_sweep.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace tpsf {
namespace {

const char* name(DecisionClassification c)
{
    switch (c) {
    case DecisionClassification::certified_satisfactory:
        return "satisfactory";
    case DecisionClassification::certified_violation:
        return "violation";
    case DecisionClassification::indeterminate:
        return "indeterminate";
    }
    return "unknown";
}

std::string optional_number(
    const std::optional<double>& value)
{
    if (!value)
        return "null";

    std::ostringstream out;
    out << std::setprecision(17) << *value;
    return out.str();
}

} // namespace

DecisionBoundarySweep assess_decision_boundary_sweep(
    const PipelineDecisionResult& baseline,
    const std::vector<double>& minimum_pressures)
{
    const auto& original = baseline.decision.points;
    const auto& time = baseline.directional.time;

    if (original.empty() || original.size() != time.size())
        throw std::invalid_argument(
            "Invalid baseline forecast history");

    DecisionBoundarySweep result;

    for (double limit : minimum_pressures) {
        if (!std::isfinite(limit))
            throw std::invalid_argument(
                "Nonfinite engineering threshold");

        BoundarySweepPoint entry;
        entry.minimum_pressure_pa = limit;

        std::vector<DecisionClassification> classes;
        classes.reserve(original.size());

        for (const auto& point : original) {
            DecisionClassification c =
                DecisionClassification::indeterminate;

            // Unbounded or invalid uncertainty must not
            // be interpreted as a finite certificate.
            if (point.finite_bounds &&
                std::isfinite(point.lower_bound) &&
                std::isfinite(point.upper_bound)) {

                if (point.lower_bound >= limit) {
                    c = DecisionClassification::
                        certified_satisfactory;
                } else if (point.upper_bound < limit) {
                    c = DecisionClassification::
                        certified_violation;
                }
            }

            classes.push_back(c);

            switch (c) {
            case DecisionClassification::certified_satisfactory:
                ++entry.satisfactory_count;
                break;
            case DecisionClassification::certified_violation:
                ++entry.violation_count;
                break;
            case DecisionClassification::indeterminate:
                ++entry.indeterminate_count;
                break;
            }
        }

        if (entry.satisfactory_count == classes.size()) {
            entry.classification =
                DecisionClassification::certified_satisfactory;
        } else if (entry.violation_count == classes.size()) {
            entry.classification =
                DecisionClassification::certified_violation;
        } else {
            entry.classification =
                DecisionClassification::indeterminate;
        }

        // Scan backward: persistent means every remaining
        // sampled time satisfies the same classification.
        bool satisfactory_suffix = true;
        bool violation_suffix = true;

        for (std::size_t k = classes.size(); k-- > 0;) {
            satisfactory_suffix =
                satisfactory_suffix &&
                classes[k] ==
                    DecisionClassification::certified_satisfactory;

            violation_suffix =
                violation_suffix &&
                classes[k] ==
                    DecisionClassification::certified_violation;

            if (satisfactory_suffix)
                entry.persistent_satisfactory_time = time[k];

            if (violation_suffix)
                entry.persistent_violation_time = time[k];
        }

        result.points.push_back(entry);
    }

    return result;
}

std::string decision_boundary_sweep_json(
    const DecisionBoundarySweep& sweep)
{
    std::ostringstream out;
    out << std::setprecision(17);

    out << "{\n";
    out << "  \"schema\": \"tpsf.decision-boundary-sweep.v1\",\n";
    out << "  \"reused_forecast\": true,\n";
    out << "  \"points\": [\n";

    for (std::size_t k = 0; k < sweep.points.size(); ++k) {
        const auto& p = sweep.points[k];

        out << "    {\n";
        out << "      \"minimum_pressure_pa\": "
            << p.minimum_pressure_pa << ",\n";
        out << "      \"classification\": \""
            << name(p.classification) << "\",\n";
        out << "      \"satisfactory_count\": "
            << p.satisfactory_count << ",\n";
        out << "      \"violation_count\": "
            << p.violation_count << ",\n";
        out << "      \"indeterminate_count\": "
            << p.indeterminate_count << ",\n";
        out << "      \"persistent_satisfactory_time_s\": "
            << optional_number(p.persistent_satisfactory_time)
            << ",\n";
        out << "      \"persistent_violation_time_s\": "
            << optional_number(p.persistent_violation_time)
            << "\n";
        out << "    }";

        if (k + 1 < sweep.points.size())
            out << ",";

        out << "\n";
    }

    out << "  ]\n}\n";
    return out.str();
}

std::string decision_boundary_sweep_markdown(
    const DecisionBoundarySweep& sweep)
{
    std::ostringstream out;
    out << std::setprecision(12);

    out << "# Decision-Boundary Sweep\n\n";
    out << "The nominal forecast and uncertainty bounds "
           "are reused for every engineering threshold.\n\n";

    out << "| Minimum pressure (MPa) | Classification | "
           "Satisfactory | Violation | Indeterminate | "
           "Persistent satisfactory (s) | "
           "Persistent violation (s) |\n";

    out << "|---:|---|---:|---:|---:|---:|---:|\n";

    for (const auto& p : sweep.points) {
        out << "| " << p.minimum_pressure_pa / 1.0e6
            << " | " << name(p.classification)
            << " | " << p.satisfactory_count
            << " | " << p.violation_count
            << " | " << p.indeterminate_count
            << " | ";

        if (p.persistent_satisfactory_time)
            out << *p.persistent_satisfactory_time;
        else
            out << "not attained";

        out << " | ";

        if (p.persistent_violation_time)
            out << *p.persistent_violation_time;
        else
            out << "not attained";

        out << " |\n";
    }

    out << "\nClassification is based on deterministic "
           "local-linearized uncertainty bounds at sampled "
           "forecast times. Persistent classifications apply "
           "only to the remaining sampled horizon. "
           "No nonlinear certificate is implied.\n";

    return out.str();
}

void write_decision_boundary_sweep(
    const DecisionBoundarySweep& sweep,
    const std::filesystem::path& directory)
{
    const auto json = decision_boundary_sweep_json(sweep);
    const auto markdown =
        decision_boundary_sweep_markdown(sweep);

    std::filesystem::create_directories(directory);

    std::ofstream json_file(
        directory / "decision_boundary_sweep.json");

    std::ofstream markdown_file(
        directory / "decision_boundary_sweep.md");

    if (!json_file || !markdown_file)
        throw std::runtime_error(
            "Cannot create decision sweep reports");

    json_file << json;
    markdown_file << markdown;

    if (!json_file || !markdown_file)
        throw std::runtime_error(
            "Cannot write decision sweep reports");
}

} // namespace tpsf
