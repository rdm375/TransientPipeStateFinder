#include "tpsf/decision_boundary_sweep.hpp"

#include <cassert>
#include <string>

int main()
{
    using namespace tpsf;

    PipelineDecisionResult baseline;

    baseline.directional.time = {0.0, 1.0, 2.0};

    for (int k = 0; k < 3; ++k) {
        LinearizedDecisionResult point;

        point.nominal = 5.0e6;
        point.lower_bound = 4.9e6;
        point.upper_bound = 5.1e6;
        point.finite_bounds = true;

        baseline.decision.points.push_back(point);
    }

    auto result = assess_decision_boundary_sweep(
        baseline,
        {4.0e6, 5.0e6, 6.0e6});

    assert(result.points.size() == 3);

    assert(
        result.points[0].classification ==
        DecisionClassification::certified_satisfactory);

    assert(
        result.points[1].classification ==
        DecisionClassification::indeterminate);

    assert(
        result.points[2].classification ==
        DecisionClassification::certified_violation);

    assert(
        result.points[0].persistent_satisfactory_time ==
        0.0);

    assert(
        result.points[2].persistent_violation_time ==
        0.0);

    // An unbounded point must never be certified.
    baseline.decision.points[1].finite_bounds = false;

    result = assess_decision_boundary_sweep(
        baseline, {4.0e6});

    assert(result.points[0].indeterminate_count == 1);
    assert(!result.points[0].persistent_satisfactory_time);

    const auto json = decision_boundary_sweep_json(result);
    const auto markdown =
        decision_boundary_sweep_markdown(result);

    assert(json.find("indeterminate") != std::string::npos);
    assert(markdown.find("not attained") != std::string::npos);

    assert(json == decision_boundary_sweep_json(result));

    return 0;
}
