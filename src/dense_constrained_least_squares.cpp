#include "tpsf/dense_constrained_least_squares.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace tpsf {
namespace {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

void require_finite(const Matrix& value)
{
    if (!value.allFinite())
        throw std::invalid_argument("Nonfinite matrix");
}

void require_finite(const Vector& value)
{
    if (!value.allFinite())
        throw std::invalid_argument("Nonfinite vector");
}

Eigen::Index numerical_rank(
    const Vector& singular_values,
    double tolerance)
{
    if (singular_values.size() == 0)
        return 0;

    const double maximum = singular_values[0];

    if (maximum == 0.0)
        return 0;

    Eigen::Index rank = 0;

    for (Eigen::Index i = 0;
         i < singular_values.size(); ++i) {
        if (singular_values[i] > tolerance*maximum)
            ++rank;
    }

    return rank;
}

double retained_minimum(
    const Vector& singular_values,
    Eigen::Index rank)
{
    return rank > 0 ? singular_values[rank-1] : 0.0;
}

} // namespace

DenseConstrainedLeastSquaresResult
solve_dense_constrained_least_squares(
    const Matrix& A,
    const Vector& c,
    const Matrix& B,
    const Vector& r,
    const DenseConstrainedLeastSquaresOptions& options)
{
    if (A.rows() != c.size() ||
        B.rows() != r.size() ||
        A.cols() != B.cols()) {
        throw std::invalid_argument(
            "Incompatible constrained least-squares dimensions");
    }

    if (!std::isfinite(options.relative_rank_tolerance) ||
        !(options.relative_rank_tolerance > 0.0) ||
        options.relative_rank_tolerance >= 1.0 ||
        !std::isfinite(options.feasibility_tolerance) ||
        !(options.feasibility_tolerance >= 0.0)) {
        throw std::invalid_argument(
            "Invalid constrained least-squares tolerances");
    }

    require_finite(A);
    require_finite(c);
    require_finite(B);
    require_finite(r);

    const Eigen::Index n = A.cols();

    DenseConstrainedLeastSquaresResult result;
    result.step = Vector::Zero(n);

    Matrix Z = Matrix::Identity(n, n);

    if (A.rows() > 0 && n > 0) {
        Eigen::JacobiSVD<Matrix> svd(
            A, Eigen::ComputeFullU | Eigen::ComputeFullV);

        const auto singular = svd.singularValues();

        const Eigen::Index rank = numerical_rank(
            singular, options.relative_rank_tolerance);

        result.constraint_rank =
            static_cast<std::size_t>(rank);

        if (singular.size() > 0)
            result.constraint_singular_value_max = singular[0];

        result.constraint_singular_value_min_retained =
            retained_minimum(singular, rank);

        if (rank > 0) {
            const Vector coefficients =
                svd.matrixU().leftCols(rank).transpose()*(-c);

            const Vector scaled =
                coefficients.array() /
                singular.head(rank).array();

            result.step =
                svd.matrixV().leftCols(rank)*scaled;
        }

        Z = svd.matrixV().rightCols(n-rank);
    }

    const double particular_defect =
        (A*result.step+c).norm();

    const double compatibility_scale =
        std::max({
            1.0,
            c.norm(),
            A.norm()*result.step.norm()
        });

    result.linearized_constraints_compatible =
        std::isfinite(particular_defect) &&
        particular_defect <=
            options.feasibility_tolerance*compatibility_scale;

    // Do not optimize an inconsistent linearization
    // as though it were an equality-constrained solution.
    if (!result.linearized_constraints_compatible) {
        result.linearized_constraint_residual =
            particular_defect;

        result.linearized_measurement_residual =
            (B*result.step+r).norm();

        return result;
    }

    const Eigen::Index nullity = Z.cols();

    if (nullity > 0 && B.rows() > 0) {
        const Matrix reduced = B*Z;
        const Vector target = -(B*result.step+r);

        Eigen::JacobiSVD<Matrix> svd(
            reduced,
            Eigen::ComputeThinU | Eigen::ComputeThinV);

        const auto singular = svd.singularValues();

        const Eigen::Index rank = numerical_rank(
            singular, options.relative_rank_tolerance);

        result.reduced_measurement_rank =
            static_cast<std::size_t>(rank);

        if (singular.size() > 0)
            result.reduced_measurement_singular_value_max =
                singular[0];

        result.reduced_measurement_singular_value_min_retained =
            retained_minimum(singular, rank);

        if (rank > 0) {
            const Vector coefficients =
                svd.matrixU().leftCols(rank).transpose()*target;

            const Vector scaled =
                coefficients.array() /
                singular.head(rank).array();

            result.step +=
                Z*(svd.matrixV().leftCols(rank)*scaled);
        }
    }

    result.unresolved_directions =
        static_cast<std::size_t>(nullity) -
        result.reduced_measurement_rank;

    result.linearized_constraint_residual =
        (A*result.step+c).norm();

    result.linearized_measurement_residual =
        (B*result.step+r).norm();

    if (!result.step.allFinite() ||
        !std::isfinite(result.linearized_constraint_residual) ||
        !std::isfinite(result.linearized_measurement_residual)) {
        throw std::runtime_error(
            "Nonfinite constrained least-squares solution");
    }

    return result;
}

} // namespace tpsf
