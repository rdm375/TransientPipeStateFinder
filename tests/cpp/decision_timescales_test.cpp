#include "tpsf/decision_timescales.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace tpsf;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

ForecastDecisionResult make_forecast(
    const std::vector<DecisionClassification>& classes)
{
    ForecastDecisionResult forecast;

    for (const auto classification : classes) {
        LinearizedDecisionResult point;
        point.classification = classification;
        forecast.points.push_back(point);
    }

    return forecast;
}

void test_persistent_satisfactory()
{
    const std::vector<double> time{0, 1, 2, 3, 4};

    const auto forecast = make_forecast({
        DecisionClassification::indeterminate,
        DecisionClassification::certified_satisfactory,
        DecisionClassification::indeterminate,
        DecisionClassification::certified_satisfactory,
        DecisionClassification::certified_satisfactory
    });

    const auto result =
        assess_decision_timescales(time, forecast);

    require(result.satisfactory.time.has_value() &&
            *result.satisfactory.time == 3,
            "Incorrect persistent satisfactory time");

    require(result.resolved.time.has_value() &&
            *result.resolved.time == 3,
            "Incorrect resolution time");

    require(result.indeterminate_count == 2,
            "Incorrect indeterminate count");

    require(result.first_indeterminate_time == 0 &&
            result.last_indeterminate_time == 2,
            "Incorrect indeterminate horizon");
}

void test_persistent_violation()
{
    const std::vector<double> time{0, 1, 2, 3};

    const auto forecast = make_forecast({
        DecisionClassification::indeterminate,
        DecisionClassification::certified_violation,
        DecisionClassification::certified_violation,
        DecisionClassification::certified_violation
    });

    const auto result =
        assess_decision_timescales(time, forecast);

    require(result.violation.time.has_value() &&
            *result.violation.time == 1,
            "Incorrect persistent violation time");

    require(!result.satisfactory.time.has_value(),
            "Unexpected satisfactory suffix");

    require(result.resolved.time.has_value() &&
            *result.resolved.time == 1,
            "Incorrect resolved suffix");
}

void test_mixed_resolved_suffix()
{
    const std::vector<double> time{0, 1, 2, 3};

    const auto forecast = make_forecast({
        DecisionClassification::indeterminate,
        DecisionClassification::certified_satisfactory,
        DecisionClassification::certified_violation,
        DecisionClassification::certified_satisfactory
    });

    const auto result =
        assess_decision_timescales(time, forecast);

    require(result.resolved.time.has_value() &&
            *result.resolved.time == 1,
            "Incorrect mixed resolved suffix");

    require(result.satisfactory.time.has_value() &&
            *result.satisfactory.time == 3,
            "Incorrect final satisfactory suffix");

    require(!result.violation.time.has_value(),
            "Unexpected persistent violation");
}

void test_unresolved_terminal()
{
    const std::vector<double> time{0, 1, 2};

    const auto forecast = make_forecast({
        DecisionClassification::certified_satisfactory,
        DecisionClassification::certified_violation,
        DecisionClassification::indeterminate
    });

    const auto result =
        assess_decision_timescales(time, forecast);

    require(!result.resolved.time.has_value(),
            "Terminal indeterminacy must prevent resolution");

    require(!result.satisfactory.time.has_value(),
            "Unexpected satisfactory suffix");

    require(!result.violation.time.has_value(),
            "Unexpected violation suffix");
}

void test_operational_adapter()
{
    PipelineDecisionResult pipeline;

    pipeline.directional.time = {0, 1, 2, 3};

    pipeline.directional.measurement_unit_radius =
        {10, 5, 1, 0.1};

    pipeline.directional.unresolved_unit_radius =
        {0, 10, 1, 0.1};

    pipeline.decision = make_forecast({
        DecisionClassification::indeterminate,
        DecisionClassification::indeterminate,
        DecisionClassification::certified_satisfactory,
        DecisionClassification::certified_satisfactory
    });

    SensitivityDecayOptions options;
    options.relative_threshold = 0.1;
    options.absolute_threshold = 0.5;

    const auto result =
        assess_operational_timescales(pipeline, options);

    require(result.decision.resolved.time == 2,
            "Incorrect operational resolution time");

    require(result.sensitivity.measurement.peak_relative.time == 2,
            "Incorrect measurement attenuation time");

    require(result.sensitivity.unresolved.peak_relative.time == 2,
            "Incorrect unresolved attenuation time");

    require(result.measurement_resolution_lag == 0,
            "Incorrect measurement resolution lag");

    require(result.unresolved_resolution_lag == 0,
            "Incorrect unresolved resolution lag");
}

void test_invalid_dimensions()
{
    const std::vector<double> time{0, 1};

    const auto forecast = make_forecast({
        DecisionClassification::indeterminate
    });

    bool threw = false;

    try {
        (void)assess_decision_timescales(time, forecast);
    } catch (const std::invalid_argument&) {
        threw = true;
    }

    require(threw, "Dimension mismatch was not rejected");
}

} // namespace

int main()
{
    try {
        test_persistent_satisfactory();
        test_persistent_violation();
        test_mixed_resolved_suffix();
        test_unresolved_terminal();
        test_operational_adapter();
        test_invalid_dimensions();

        std::cout
            << "PASS E5.6b operational decision timescales\n";

        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL E5.6b: "
                  << error.what() << '\n';
        return 1;
    }
}
