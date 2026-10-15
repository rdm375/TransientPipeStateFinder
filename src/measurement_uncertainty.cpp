#include "tpsf/measurement_uncertainty.hpp"

#include "tpsf/dense_full_space_subproblem.hpp"

#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace tpsf {
namespace {

void validate_options(const MeasurementUncertaintyOptions& options)
{
    if (!std::isfinite(options.relative_rank_tolerance) ||
        options.relative_rank_tolerance <= 0.0 ||
        options.relative_rank_tolerance >= 1.0 ||
        !std::isfinite(options.unresolved_relative_tolerance) ||
        options.unresolved_relative_tolerance < 0.0 ||
        !std::isfinite(options.unresolved_absolute_tolerance) ||
        options.unresolved_absolute_tolerance < 0.0)
        throw std::invalid_argument(
            "Invalid measurement uncertainty options");
}

} // namespace

EngineeringUncertainty propagate_measurement_uncertainty(
    const Eigen::MatrixXd& S,
    const Eigen::VectorXd& g,
    const MeasurementUncertaintyOptions& options)
{
    validate_options(options);

    if (S.cols() != g.size() ||
        !S.allFinite() ||
        !g.allFinite())
        throw std::invalid_argument(
            "Invalid uncertainty matrix or gradient");

    EngineeringUncertainty result;

    const Eigen::Index n = S.cols();

    result.feasible_dimension =
        static_cast<std::size_t>(n);

    // No feasible degrees of freedom: the quantity is locally
    // fixed by the equality constraints.
    if (n == 0) {
        result.status = UncertaintyStatus::finite;
        result.finite_data_only_bound = true;
        return result;
    }

    // No observations: every feasible direction is unresolved.
    if (S.rows() == 0) {
        result.unresolved_directions =
            static_cast<std::size_t>(n);

        result.unresolved_sensitivity_norm = g.norm();

        const double threshold =
            options.unresolved_absolute_tolerance +
            options.unresolved_relative_tolerance * g.norm();

        if (g.norm() <= threshold) {
            result.status = UncertaintyStatus::finite;
            result.finite_data_only_bound = true;
        } else {
            result.status = UncertaintyStatus::unbounded;
            result.finite_data_only_bound = false;
            result.variance =
                std::numeric_limits<double>::infinity();
            result.standard_deviation =
                std::numeric_limits<double>::infinity();
        }

        return result;
    }

    // Full V is necessary when S is rectangular and
    // rank deficient.
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        S, Eigen::ComputeThinU | Eigen::ComputeFullV);

    const auto singular = svd.singularValues();

    Eigen::Index rank = 0;

    const double largest =
        singular.size() > 0 ? singular[0] : 0.0;

    if (largest > 0.0) {
        for (Eigen::Index i = 0;
             i < singular.size(); ++i) {
            if (singular[i] >
                options.relative_rank_tolerance * largest)
                ++rank;
        }
    }

    result.information_rank =
        static_cast<std::size_t>(rank);

    result.unresolved_directions =
        static_cast<std::size_t>(n-rank);

    const Eigen::MatrixXd V = svd.matrixV();

    const Eigen::VectorXd components =
        V.transpose() * g;

    const double unresolved_norm =
        components.tail(n-rank).norm();

    result.unresolved_sensitivity_norm =
        unresolved_norm;

    const double threshold =
        options.unresolved_absolute_tolerance +
        options.unresolved_relative_tolerance * g.norm();

    if (unresolved_norm > threshold) {
        result.status = UncertaintyStatus::unbounded;
        result.finite_data_only_bound = false;
        result.variance =
            std::numeric_limits<double>::infinity();
        result.standard_deviation =
            std::numeric_limits<double>::infinity();

        return result;
    }

    // Only retained singular directions contribute.
    // Standardized measurement errors have unit variance.
    double variance = 0.0;

    for (Eigen::Index i = 0; i < rank; ++i) {
        const double contribution =
            components[i] / singular[i];

        variance += contribution * contribution;
    }

    if (!std::isfinite(variance))
        throw std::runtime_error(
            "Nonfinite propagated measurement variance");

    result.variance = variance;
    result.standard_deviation = std::sqrt(variance);
    result.status = UncertaintyStatus::finite;
    result.finite_data_only_bound = true;

    return result;
}

EngineeringUncertainty assess_engineering_uncertainty(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& coordinates,
    const StatefindingQuality& quality,
    const QuantitySensitivity& sensitivity,
    const MeasurementUncertaintyOptions& options)
{
    if (!sensitivity.differentiable) {
        EngineeringUncertainty result;
        result.status = UncertaintyStatus::nondifferentiable;
        return result;
    }

    const auto assembled =
        assemble_dense_full_space_subproblem(
            problem, coordinates);

    const Eigen::Index n =
        static_cast<Eigen::Index>(
            problem.free_variable_count());

    if (quality.feasible_basis.rows() != n ||
        sensitivity.gradient.size() != n)
        throw std::invalid_argument(
            "Incompatible quality and sensitivity dimensions");

    const Eigen::MatrixXd S =
        assembled.measurement_jacobian *
        quality.feasible_basis;

    const Eigen::VectorXd g =
        quality.feasible_basis.transpose() *
        sensitivity.gradient;

    return propagate_measurement_uncertainty(
        S, g, options);
}

} // namespace tpsf
