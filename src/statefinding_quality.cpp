#include "tpsf/statefinding_quality.hpp"
#include "tpsf/dense_full_space_subproblem.hpp"

#include <Eigen/SVD>

#include <cmath>
#include <stdexcept>

namespace tpsf {
namespace {

Eigen::Index numerical_rank(
    const Eigen::VectorXd& singular,
    double tolerance)
{
    if (singular.size() == 0 || singular[0] <= 0.0)
        return 0;

    Eigen::Index rank = 0;

    for (Eigen::Index i = 0; i < singular.size(); ++i) {
        if (singular[i] > tolerance * singular[0])
            ++rank;
    }

    return rank;
}

} // namespace

StatefindingQuality assess_statefinding_quality(
    const FullSpaceProblem& problem,
    const DenseSQPResult& solution,
    const StatefindingQualityOptions& options)
{
    if (!std::isfinite(options.relative_rank_tolerance) ||
        options.relative_rank_tolerance <= 0.0 ||
        options.relative_rank_tolerance >= 1.0)
        throw std::invalid_argument(
            "Invalid relative rank tolerance");

    const auto assembled =
        assemble_dense_full_space_subproblem(
            problem, solution.coordinates);

    const auto& A = assembled.constraint_jacobian;
    const auto& B = assembled.measurement_jacobian;

    const Eigen::Index n = A.cols();

    StatefindingQuality result;

    result.termination = solution.termination;
    result.converged = solution.converged();

    result.objective =
        0.5 * assembled.measurement_residual.squaredNorm();

    result.constraint_l2 =
        assembled.constraint_residual.norm();

    result.measurement_l2 =
        assembled.measurement_residual.norm();

    result.variables = static_cast<std::size_t>(n);
    result.constraints =
        static_cast<std::size_t>(A.rows());
    result.measurements =
        static_cast<std::size_t>(B.rows());

    const Eigen::VectorXd gradient =
        B.transpose() * assembled.measurement_residual;

    // Full V is essential: thin V omits nullspace directions
    // for rectangular constraint Jacobians.
    Eigen::JacobiSVD<Eigen::MatrixXd> constraint_svd(
        A, Eigen::ComputeFullV);

    const Eigen::Index rank_a =
        numerical_rank(
            constraint_svd.singularValues(),
            options.relative_rank_tolerance);

    result.constraint_rank =
        static_cast<std::size_t>(rank_a);

    const Eigen::Index feasible_dimension = n - rank_a;

    result.feasible_dimension =
        static_cast<std::size_t>(feasible_dimension);

    result.feasible_basis =
        constraint_svd.matrixV().rightCols(
            feasible_dimension);

    // At a feasible point, Z^T gradient measures
    // first-order stationarity along feasible directions.
    result.projected_gradient =
        (result.feasible_basis.transpose() * gradient).norm();

    if (feasible_dimension == 0) {
        result.unresolved_basis.resize(n, 0);
        result.information_singular_values.resize(0);
        result.locally_identifiable = true;
        return result;
    }

    const Eigen::MatrixXd S =
        B * result.feasible_basis;

    // Explicitly handle zero-measurement cases: Eigen SVD
    // does not support every empty-matrix configuration.
    if (S.rows() == 0) {
        result.information_singular_values =
            Eigen::VectorXd::Zero(feasible_dimension);

        result.unresolved_basis = result.feasible_basis;
        result.unresolved_directions =
            static_cast<std::size_t>(feasible_dimension);

        return result;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> information_svd(
        S, Eigen::ComputeFullV);

    const Eigen::Index rank_s =
        numerical_rank(
            information_svd.singularValues(),
            options.relative_rank_tolerance);

    const Eigen::Index unresolved =
        feasible_dimension - rank_s;

    result.measurement_rank =
        static_cast<std::size_t>(rank_s);

    result.unresolved_directions =
        static_cast<std::size_t>(unresolved);

    result.information_singular_values =
        Eigen::VectorXd::Zero(feasible_dimension);

    result.information_singular_values.head(
        information_svd.singularValues().size()) =
            information_svd.singularValues();

    result.unresolved_basis =
        result.feasible_basis *
        information_svd.matrixV().rightCols(unresolved);

    result.locally_identifiable = (unresolved == 0);

    return result;
}

} // namespace tpsf
