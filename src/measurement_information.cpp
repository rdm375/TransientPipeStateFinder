#include "tpsf/measurement_information.hpp"

#include <Eigen/SVD>

#include <cmath>
#include <stdexcept>

namespace tpsf {

MeasurementInformationDecomposition
decompose_measurement_information(
    const Eigen::MatrixXd& B,
    const Eigen::MatrixXd& Z,
    double relative_rank_tolerance)
{
    if (!std::isfinite(relative_rank_tolerance) ||
        relative_rank_tolerance <= 0.0 ||
        relative_rank_tolerance >= 1.0)
        throw std::invalid_argument(
            "Invalid information rank tolerance");

    if (B.cols() != Z.rows() ||
        !B.allFinite() ||
        !Z.allFinite())
        throw std::invalid_argument(
            "Invalid information decomposition dimensions");

    const Eigen::Index n = Z.rows();
    const Eigen::Index d = Z.cols();

    MeasurementInformationDecomposition result;
    result.feasible_dimension =
        static_cast<std::size_t>(d);

    result.resolved_basis.resize(n, 0);
    result.unresolved_basis.resize(n, 0);
    result.retained_singular_values.resize(0);

    if (d == 0)
        return result;

    // A feasible basis must be orthonormal.
    const Eigen::MatrixXd gram =
        Z.transpose() * Z;

    const Eigen::MatrixXd identity =
        Eigen::MatrixXd::Identity(d, d);

    if ((gram - identity).norm() >
        1.0e-10 * static_cast<double>(d))
        throw std::invalid_argument(
            "Feasible basis is not orthonormal");

    if (B.rows() == 0) {
        result.unresolved_basis = Z;
        result.unresolved_dimension =
            static_cast<std::size_t>(d);
        return result;
    }

    const Eigen::MatrixXd S = B * Z;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        S, Eigen::ComputeFullV);

    const Eigen::VectorXd singular =
        svd.singularValues();

    const double largest =
        singular.size() > 0 ? singular[0] : 0.0;

    Eigen::Index rank = 0;

    if (largest > 0.0) {
        for (Eigen::Index j = 0;
             j < singular.size(); ++j) {
            if (singular[j] >
                relative_rank_tolerance * largest)
                ++rank;
        }
    }

    const Eigen::Index unresolved = d - rank;

    const Eigen::MatrixXd V = svd.matrixV();

    result.resolved_basis =
        Z * V.leftCols(rank);

    result.unresolved_basis =
        Z * V.rightCols(unresolved);

    result.retained_singular_values =
        singular.head(rank);

    result.information_rank =
        static_cast<std::size_t>(rank);

    result.unresolved_dimension =
        static_cast<std::size_t>(unresolved);

    return result;
}

} // namespace tpsf
