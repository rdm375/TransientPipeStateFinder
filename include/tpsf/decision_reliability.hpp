#pragma once

#include <Eigen/Dense>

#include <cstddef>
#include <limits>
#include <vector>

namespace tpsf {

enum class DecisionLimit {
    minimum,
    maximum
};

enum class DecisionClassification {
    certified_satisfactory,
    certified_violation,
    indeterminate
};

// A local, linearized uncertainty model:
//
//     delta_u = L*z + delta_unbounded, ||z||_2 <= 1.
//
// The columns of L represent bounded, possibly correlated
// perturbations in the same coordinates as the gradient.
//
// The caller must explicitly identify whether the quantity
// depends on any unbounded uncertainty directions.
struct LinearizedDecisionInput {
    double nominal = 0.0;
    double limit = 0.0;

    DecisionLimit limit_type = DecisionLimit::minimum;

    Eigen::VectorXd gradient;
    Eigen::MatrixXd uncertainty_factor;

    bool sensitive_to_unbounded_directions = false;
};

struct LinearizedDecisionResult {
    DecisionClassification classification =
        DecisionClassification::indeterminate;

    double nominal = 0.0;
    double limit = 0.0;

    double lower_bound =
        -std::numeric_limits<double>::infinity();

    double upper_bound =
        std::numeric_limits<double>::infinity();

    double uncertainty_radius =
        std::numeric_limits<double>::infinity();

    // Positive means the entire uncertainty interval
    // satisfies the engineering requirement.
    double worst_case_margin =
        -std::numeric_limits<double>::infinity();

    bool finite_bounds = false;
};

LinearizedDecisionResult assess_linearized_decision(
    const LinearizedDecisionInput& input);

struct ForecastDecisionResult {
    std::vector<LinearizedDecisionResult> points;

    // The entire assessed interval is satisfactory only
    // when every point is certified satisfactory.
    DecisionClassification classification =
        DecisionClassification::indeterminate;

    // First point after which all remaining points are
    // certified satisfactory. The index is absent if no
    // such suffix exists.
    bool has_reliable_suffix = false;
    std::size_t reliable_suffix_start = 0;
};

ForecastDecisionResult assess_forecast_decisions(
    const std::vector<LinearizedDecisionInput>& inputs);

} // namespace tpsf
