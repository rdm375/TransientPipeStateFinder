#include "tpsf/decision_reliability.hpp"

#include <Eigen/Dense>

#include <cmath>
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

void near(double actual, double expected, double tol,
          const std::string& message)
{
    require(std::abs(actual-expected) <= tol, message);
}

LinearizedDecisionInput example(
    double nominal,
    double limit,
    DecisionLimit type)
{
    LinearizedDecisionInput input;
    input.nominal = nominal;
    input.limit = limit;
    input.limit_type = type;

    input.gradient = Eigen::Vector2d(3.0, 4.0);
    input.uncertainty_factor =
        Eigen::Matrix2d::Identity();

    return input;
}

void test_ellipsoid()
{
    auto input = example(
        100.0, 94.0, DecisionLimit::minimum);

    const auto result = assess_linearized_decision(input);

    near(result.uncertainty_radius, 5.0, 1e-12,
         "Incorrect ellipsoid support function");
    near(result.lower_bound, 95.0, 1e-12,
         "Incorrect lower bound");
    near(result.upper_bound, 105.0, 1e-12,
         "Incorrect upper bound");

    require(result.classification ==
                DecisionClassification::certified_satisfactory,
            "Expected satisfactory lower limit");
}

void test_correlated_uncertainty()
{
    auto input = example(
        100.0, 90.0, DecisionLimit::minimum);

    input.uncertainty_factor <<
        2.0, 1.0,
        0.0, 3.0;

    // L^T*g = (6,15).
    const auto result = assess_linearized_decision(input);

    near(result.uncertainty_radius,
         std::sqrt(261.0), 1e-12,
         "Incorrect correlated uncertainty radius");
}

void test_three_way_classification()
{
    auto input = example(
        100.0, 106.0, DecisionLimit::minimum);

    require(assess_linearized_decision(input).classification ==
                DecisionClassification::certified_violation,
            "Expected certified lower-limit violation");

    input.limit = 102.0;

    require(assess_linearized_decision(input).classification ==
                DecisionClassification::indeterminate,
            "Expected indeterminate lower-limit decision");

    input.limit_type = DecisionLimit::maximum;
    input.limit = 106.0;

    require(assess_linearized_decision(input).classification ==
                DecisionClassification::certified_satisfactory,
            "Expected satisfactory upper limit");

    input.limit = 94.0;

    require(assess_linearized_decision(input).classification ==
                DecisionClassification::certified_violation,
            "Expected certified upper-limit violation");

    input.limit = 102.0;

    require(assess_linearized_decision(input).classification ==
                DecisionClassification::indeterminate,
            "Expected indeterminate upper-limit decision");
}

void test_unbounded()
{
    auto input = example(
        100.0, 1.0, DecisionLimit::minimum);

    input.sensitive_to_unbounded_directions = true;

    const auto result = assess_linearized_decision(input);

    require(!result.finite_bounds,
            "Unbounded uncertainty must not yield finite bounds");

    require(result.classification ==
                DecisionClassification::indeterminate,
            "Unbounded uncertainty must be indeterminate");
}

// E5.5d: strict unbounded-sensitivity regression
void test_tiny_unbounded_sensitivity()
{
    // The unresolved coordinate is unbounded.
    //
    // Even an extremely small nonzero derivative
    // can therefore produce an arbitrarily large
    // change in the local linearized forecast.

    const double sensitivity = 1.0e-14;

    auto input = example(
        100.0, 90.0, DecisionLimit::minimum);

    // Isolate the unbounded-direction classification
    // from bounded measurement uncertainty.
    input.gradient = Eigen::VectorXd::Zero(1);
    input.uncertainty_factor =
        Eigen::MatrixXd::Zero(1, 1);

    input.sensitive_to_unbounded_directions =
        sensitivity > 0.0;

    const auto result =
        assess_linearized_decision(input);

    require(
        !result.finite_bounds,
        "Tiny nonzero unbounded sensitivity must not "
        "produce finite bounds");

    require(
        result.classification ==
            DecisionClassification::indeterminate,
        "Tiny nonzero unbounded sensitivity must "
        "remain indeterminate");

    // Exactly zero sensitivity is different.
    input.sensitive_to_unbounded_directions = false;

    const auto zero_result =
        assess_linearized_decision(input);

    require(
        zero_result.finite_bounds,
        "Zero unresolved sensitivity should permit "
        "finite bounds");

    require(
        zero_result.classification ==
            DecisionClassification::certified_satisfactory,
        "Zero unresolved sensitivity should permit "
        "a satisfactory classification");

    std::cout
        << "PASS tiny unbounded sensitivity regression\n";
}

void test_forecast_suffix()
{
    std::vector<LinearizedDecisionInput> inputs;

    // Unreliable at the first sample, then reliably
    // satisfactory for the remainder of the forecast.
    inputs.push_back(example(
        100.0, 102.0, DecisionLimit::minimum));

    inputs.push_back(example(
        100.0, 94.0, DecisionLimit::minimum));

    inputs.push_back(example(
        100.0, 94.0, DecisionLimit::minimum));

    const auto result = assess_forecast_decisions(inputs);

    require(result.classification ==
                DecisionClassification::indeterminate,
            "Mixed interval must be indeterminate");

    require(result.has_reliable_suffix,
            "Expected reliable forecast suffix");

    require(result.reliable_suffix_start == 1,
            "Incorrect reliable suffix start");
}


void test_single_interval_violation()
{
    using namespace tpsf;

    std::vector<LinearizedDecisionInput> inputs{
        example(100.0, 94.0, DecisionLimit::minimum),
        example(100.0, 106.0, DecisionLimit::minimum),
        example(100.0, 94.0, DecisionLimit::minimum)
    };

    const auto result = assess_forecast_decisions(inputs);

    require(
        result.classification ==
            DecisionClassification::certified_violation,
        "One certified violation must invalidate an all-times requirement");

    require(result.has_reliable_suffix,
            "Expected satisfactory final suffix");

    require(result.reliable_suffix_start == 2,
            "Incorrect suffix after violation");
}

void test_invalid_inputs()
{
    auto input = example(
        100.0, 90.0, DecisionLimit::minimum);

    input.uncertainty_factor =
        Eigen::MatrixXd::Identity(3, 3);

    bool threw = false;

    try {
        (void)assess_linearized_decision(input);
    } catch (const std::invalid_argument&) {
        threw = true;
    }

    require(threw, "Expected dimension validation");
}

} // namespace

int main()
{
    try {
        test_ellipsoid();
        test_correlated_uncertainty();
        test_three_way_classification();
        test_unbounded();
        test_tiny_unbounded_sensitivity();
        test_forecast_suffix();
        test_single_interval_violation();
        test_invalid_inputs();

        std::cout
            << "PASS E5.5a linearized engineering decisions\n";

        return 0;
    } catch (const std::exception& error) {
        std::cerr
            << "FAIL E5.5a: " << error.what() << '\n';

        return 1;
    }
}
