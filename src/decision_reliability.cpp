#include "tpsf/decision_reliability.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tpsf {
namespace {

void validate(const LinearizedDecisionInput& input)
{
    if (!std::isfinite(input.nominal) ||
        !std::isfinite(input.limit))
        throw std::invalid_argument(
            "Nonfinite engineering quantity or limit");

    if (input.gradient.size() < 1 ||
        !input.gradient.allFinite())
        throw std::invalid_argument(
            "Invalid engineering gradient");

    if (input.uncertainty_factor.rows() !=
            input.gradient.size() ||
        !input.uncertainty_factor.allFinite())
        throw std::invalid_argument(
            "Invalid uncertainty factor");

    if (input.limit_type != DecisionLimit::minimum &&
        input.limit_type != DecisionLimit::maximum)
        throw std::invalid_argument(
            "Invalid engineering limit type");
}

} // namespace

LinearizedDecisionResult assess_linearized_decision(
    const LinearizedDecisionInput& input)
{
    validate(input);

    LinearizedDecisionResult result;
    result.nominal = input.nominal;
    result.limit = input.limit;

    if (!input.sensitive_to_unbounded_directions) {
        // Exact support function of the ellipsoid:
        //
        // max_{||z||<=1} g^T L z = ||L^T g||_2.
        //
        // A factor with zero columns represents a
        // deterministic quantity in the specified model.
        const double radius =
            (input.uncertainty_factor.transpose() *
             input.gradient).norm();

        if (!std::isfinite(radius))
            throw std::overflow_error(
                "Nonfinite uncertainty radius");

        result.uncertainty_radius = radius;
        result.lower_bound = input.nominal-radius;
        result.upper_bound = input.nominal+radius;

        if (!std::isfinite(result.lower_bound) ||
            !std::isfinite(result.upper_bound))
            throw std::overflow_error(
                "Nonfinite decision bounds");

        result.finite_bounds = true;
    }

    if (input.limit_type == DecisionLimit::minimum) {
        result.worst_case_margin =
            result.lower_bound-input.limit;

        if (result.lower_bound >= input.limit)
            result.classification =
                DecisionClassification::certified_satisfactory;
        else if (result.upper_bound < input.limit)
            result.classification =
                DecisionClassification::certified_violation;
    } else {
        result.worst_case_margin =
            input.limit-result.upper_bound;

        if (result.upper_bound <= input.limit)
            result.classification =
                DecisionClassification::certified_satisfactory;
        else if (result.lower_bound > input.limit)
            result.classification =
                DecisionClassification::certified_violation;
    }

    return result;
}

ForecastDecisionResult assess_forecast_decisions(
    const std::vector<LinearizedDecisionInput>& inputs)
{
    if (inputs.empty())
        throw std::invalid_argument(
            "Cannot assess an empty forecast");

    ForecastDecisionResult result;
    result.points.reserve(inputs.size());

    bool all_satisfactory = true;
    bool any_violation = false;

    for (const auto& input : inputs) {
        const auto point = assess_linearized_decision(input);
        result.points.push_back(point);

        all_satisfactory &=
            point.classification ==
            DecisionClassification::certified_satisfactory;

        any_violation |=
            point.classification ==
            DecisionClassification::certified_violation;
    }

    if (all_satisfactory)
        result.classification =
            DecisionClassification::certified_satisfactory;
    else if (any_violation)
        result.classification =
            DecisionClassification::certified_violation;

    // A single certified violation disproves an
    // all-times requirement. Otherwise a mixed interval
    // remains indeterminate.
    //
    // Search backwards for a continuously satisfactory
    // suffix of the sampled forecast.
    bool suffix_satisfactory = true;

    for (std::size_t i = result.points.size(); i-- > 0;) {
        suffix_satisfactory &=
            result.points[i].classification ==
            DecisionClassification::certified_satisfactory;

        if (suffix_satisfactory) {
            result.has_reliable_suffix = true;
            result.reliable_suffix_start = i;
        }
    }

    return result;
}

} // namespace tpsf
