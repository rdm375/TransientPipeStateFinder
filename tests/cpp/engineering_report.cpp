#include "tpsf/engineering_report.hpp"

#include <cassert>
#include <stdexcept>
#include <string>

int main()
{
    using namespace tpsf;

    EngineeringReport report;

    report.scenario.cells = 2;
    report.scenario.reconstruction_steps = 2;
    report.scenario.forecast_steps = 1;
    report.scenario.dt = 1.0;
    report.scenario.theta = 0.65;
    report.scenario.minimum_outlet_pressure_pa = 4.0e6;

    report.quality.feasible_dimension = 2;
    report.quality.measurement_rank = 1;
    report.quality.unresolved_directions = 1;
    report.quality.information_singular_values =
        Eigen::Vector2d(2.0, 0.0);

    report.decision.decision.classification =
        DecisionClassification::indeterminate;

    LinearizedDecisionResult point;
    point.classification =
        DecisionClassification::indeterminate;
    point.nominal = 5.0e6;
    point.limit = 4.0e6;
    point.finite_bounds = false;

    report.decision.decision.points.push_back(point);
    report.decision.directional.time.push_back(0.0);
    report.decision.directional.nominal_outlet_pressure
        .push_back(5.0e6);

    const auto json = engineering_report_json(report);
    const auto markdown =
        engineering_report_markdown(report);

    assert(json.find(
        "\"classification\": \"indeterminate\"")
        != std::string::npos);

    assert(json.find(
        "\"information_condition_finite\": false")
        != std::string::npos);

    assert(json.find(
        "\"statefinder_revision\": \"unknown\"")
        != std::string::npos);

    assert(json.find(
        "\"lower_bound_pa\": null")
        != std::string::npos);

    assert(markdown.find(
        "Unresolved directions: 1")
        != std::string::npos);

    assert(markdown.find(
        "not a statistical confidence level")
        != std::string::npos);

    // Rendering must be deterministic.
    assert(json == engineering_report_json(report));
    assert(markdown ==
           engineering_report_markdown(report));

    // Invalid available measurement uncertainty
    // must be rejected.
    EngineeringMeasurement measurement;
    measurement.available[0] = true;
    measurement.sigma[0] = 0.0;

    report.scenario.measurements.push_back(measurement);

    bool rejected = false;

    try {
        (void)engineering_report_json(report);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }

    assert(rejected);

    return 0;
}
