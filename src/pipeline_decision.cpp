#include "tpsf/pipeline_decision.hpp"

#include "tpsf/dense_full_space_subproblem.hpp"
#include "tpsf/measurement_information.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace tpsf {

PipelineDecisionResult assess_pipeline_decision(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& coordinates,
    const StatefindingQuality& quality,
    const Configuration& forecast_configuration,
    const BoundarySchedule& future_boundary,
    const PipelineDecisionOptions& options)
{
    if (!std::isfinite(options.minimum_outlet_pressure_pa) ||
        !std::isfinite(options.measurement_radius) ||
        options.measurement_radius < 0.0 ||
        !std::isfinite(options.unresolved_radius) ||
        options.unresolved_radius < 0.0)
        throw std::invalid_argument(
            "Invalid pipeline decision options");

    if (!quality.converged)
        throw std::invalid_argument(
            "Pipeline decision requires converged reconstruction");

    const auto assembled =
        assemble_dense_full_space_subproblem(
            problem, coordinates);

    const auto information =
        decompose_measurement_information(
            assembled.measurement_jacobian,
            quality.feasible_basis,
            quality.relative_rank_tolerance);

    if (information.information_rank !=
            quality.measurement_rank ||
        information.unresolved_dimension !=
            quality.unresolved_directions)
        throw std::runtime_error(
            "Measurement decomposition rank mismatch");

    PipelineDecisionResult result;

    // Preserve E5.4 diagnostics for comparison.
    result.attenuation = assess_forecast_attenuation(
        problem,
        coordinates,
        quality,
        forecast_configuration,
        future_boundary);

    result.directional =
        assess_forecast_directional_uncertainty(
            problem,
            coordinates,
            information,
            forecast_configuration,
            future_boundary);

    result.unresolved_bounded = options.bound_unresolved;

    // Retain this legacy diagnostic, but do not use it
    // to calculate the decision uncertainty.
    if (information.information_rank > 0) {
        const double smallest =
            information.retained_singular_values[
                information.retained_singular_values.size() - 1];

        result.measurement_coordinate_radius =
            options.measurement_radius / smallest;
    }

    std::vector<LinearizedDecisionInput> inputs;
    inputs.reserve(result.directional.time.size());

    for (std::size_t k = 0;
         k < result.directional.time.size(); ++k) {

        const double measured =
            options.measurement_radius *
            result.directional.measurement_unit_radius[k];

        const double unresolved =
            options.unresolved_radius *
            result.directional.unresolved_unit_radius[k];

        const double radius =
            measured +
            (options.bound_unresolved ? unresolved : 0.0);

        if (!std::isfinite(radius))
            throw std::overflow_error(
                "Nonfinite combined decision radius");

        LinearizedDecisionInput input;

        input.nominal =
            result.directional.nominal_outlet_pressure[k];

        input.limit =
            options.minimum_outlet_pressure_pa;

        input.limit_type = DecisionLimit::minimum;

        // Exact support radius for the two independent
        // bounded uncertainty sets, represented as one
        // equivalent scalar interval.
        input.gradient = Eigen::VectorXd::Ones(1);

        input.uncertainty_factor =
            Eigen::MatrixXd::Constant(1, 1, radius);

        input.sensitive_to_unbounded_directions =
            !options.bound_unresolved &&
            result.directional.unresolved_influence[k];

        inputs.push_back(std::move(input));
    }

    result.decision = assess_forecast_decisions(inputs);

    return result;
}

} // namespace tpsf
