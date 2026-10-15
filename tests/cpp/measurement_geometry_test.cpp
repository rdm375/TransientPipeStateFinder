#include "tpsf/measurement_information.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected,
          double tolerance, const std::string& message)
{
    require(
        std::isfinite(actual) &&
        std::isfinite(expected) &&
        std::abs(actual - expected) <= tolerance,
        message
    );
}

double directional_radius(
    const tpsf::MeasurementInformationDecomposition& info,
    const Eigen::VectorXd& gradient)
{
    require(
        info.unresolved_dimension == 0,
        "Expected fully resolved feasible subspace"
    );

    double squared = 0.0;

    for (Eigen::Index j = 0;
         j < info.resolved_basis.cols(); ++j) {

        const double sensitivity =
            gradient.dot(info.resolved_basis.col(j));

        const double singular =
            info.retained_singular_values[j];

        require(singular > 0.0, "Invalid singular value");

        const double contribution = sensitivity / singular;

        squared += contribution * contribution;
    }

    return std::sqrt(squared);
}

void test_geometry()
{
    Eigen::Matrix3d B = Eigen::Matrix3d::Zero();

    B(0, 0) = 100.0;
    B(1, 1) = 10.0;
    B(2, 2) = 1.0;

    const Eigen::Matrix3d Z =
        Eigen::Matrix3d::Identity();

    Eigen::Vector3d gradient;
    gradient << 100.0, 10.0, 1.0;

    const auto info =
        tpsf::decompose_measurement_information(
            B, Z, 1.0e-12);

    require(info.information_rank == 3,
            "Incorrect information rank");

    require(info.unresolved_dimension == 0,
            "Unexpected unresolved direction");

    const double new_radius =
        directional_radius(info, gradient);

    const double smallest =
        info.retained_singular_values[
            info.retained_singular_values.size() - 1];

    const double old_radius =
        gradient.norm() / smallest;

    near(
        new_radius,
        std::sqrt(3.0),
        1.0e-10,
        "Incorrect direction-specific radius"
    );

    near(
        old_radius,
        std::sqrt(10101.0),
        1.0e-10,
        "Incorrect conservative radius"
    );

    require(
        new_radius < old_radius,
        "Directional weighting provided no improvement"
    );

    // Rotate the feasible basis and verify invariance.
    const double angle = 0.37;

    Eigen::Matrix3d rotation =
        Eigen::Matrix3d::Identity();

    rotation(0, 0) = std::cos(angle);
    rotation(0, 1) = -std::sin(angle);
    rotation(1, 0) = std::sin(angle);
    rotation(1, 1) = std::cos(angle);

    const auto rotated =
        tpsf::decompose_measurement_information(
            B, rotation, 1.0e-12);

    const double rotated_radius =
        directional_radius(rotated, gradient);

    near(
        rotated_radius,
        new_radius,
        1.0e-10,
        "Uncertainty is not basis invariant"
    );

    // Independent inverse-information calculation.
    const Eigen::Vector3d x =
        (B.transpose() * B).ldlt().solve(gradient);

    const double independent_radius =
        std::sqrt(gradient.dot(x));

    near(
        new_radius,
        independent_radius,
        1.0e-10,
        "Disagreement with inverse-information oracle"
    );

    std::cout
        << "Old radius: " << old_radius << '\n'
        << "New radius: " << new_radius << '\n'
        << "Improvement factor: "
        << old_radius / new_radius << '\n'
        << "Rotated radius: " << rotated_radius << '\n'
        << "Independent radius: "
        << independent_radius << '\n'
        << "PASS multidirectional measurement geometry\n";
}

} // namespace

int main()
{
    try {
        test_geometry();
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL measurement geometry: "
            << error.what() << '\n';
        return 1;
    }
}
