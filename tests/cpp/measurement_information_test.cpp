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

void test_full_rank()
{
    using namespace tpsf;

    Eigen::MatrixXd Z =
        Eigen::MatrixXd::Identity(3, 2);

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2, 3);
    B(0, 0) = 4.0;
    B(1, 1) = 2.0;

    const auto result =
        decompose_measurement_information(B, Z);

    require(result.information_rank == 2,
            "Incorrect full information rank");

    require(result.unresolved_dimension == 0,
            "Unexpected unresolved directions");

    require(
        std::abs(result.retained_singular_values[0] - 4.0)
            < 1.0e-12,
        "Incorrect first singular value");

    require(
        std::abs(result.retained_singular_values[1] - 2.0)
            < 1.0e-12,
        "Incorrect second singular value");
}

void test_rank_deficient()
{
    using namespace tpsf;

    Eigen::MatrixXd Z =
        Eigen::MatrixXd::Identity(3, 3);

    Eigen::MatrixXd B(2, 3);
    B << 1.0, 2.0, 0.0,
         2.0, 4.0, 0.0;

    const auto result =
        decompose_measurement_information(B, Z);

    require(result.information_rank == 1,
            "Incorrect deficient information rank");

    require(result.unresolved_dimension == 2,
            "Incorrect unresolved dimension");

    const auto& R = result.resolved_basis;
    const auto& N = result.unresolved_basis;

    require(
        (R.transpose() * N).norm() < 1.0e-12,
        "Resolved and unresolved bases not orthogonal");

    require(
        (N.transpose() * N -
         Eigen::MatrixXd::Identity(2, 2)).norm()
            < 1.0e-12,
        "Unresolved basis not orthonormal");

    const Eigen::MatrixXd projector =
        R * R.transpose() + N * N.transpose();

    require(
        (projector -
         Eigen::MatrixXd::Identity(3, 3)).norm()
            < 1.0e-12,
        "Decomposition does not span feasible space");

    require(
        (B * N).norm() < 1.0e-12,
        "Unresolved directions affect measurements");
}

void test_no_measurements()
{
    using namespace tpsf;

    Eigen::MatrixXd Z =
        Eigen::MatrixXd::Identity(3, 2);

    Eigen::MatrixXd B(0, 3);

    const auto result =
        decompose_measurement_information(B, Z);

    require(result.information_rank == 0,
            "Unexpected information rank");

    require(result.unresolved_dimension == 2,
            "Missing unresolved directions");

    require(
        (result.unresolved_basis - Z).norm() == 0.0,
        "Zero-measurement basis changed");
}

void test_no_feasible_directions()
{
    using namespace tpsf;

    Eigen::MatrixXd Z(3, 0);
    Eigen::MatrixXd B(2, 3);
    B.setZero();

    const auto result =
        decompose_measurement_information(B, Z);

    require(result.information_rank == 0,
            "Unexpected rank");

    require(result.unresolved_dimension == 0,
            "Unexpected unresolved dimension");
}

void test_invalid_basis()
{
    using namespace tpsf;

    Eigen::MatrixXd Z(2, 1);
    Z << 2.0, 0.0;

    Eigen::MatrixXd B =
        Eigen::MatrixXd::Identity(2, 2);

    bool rejected = false;

    try {
        (void)decompose_measurement_information(B, Z);
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    require(rejected,
            "Nonorthonormal basis was accepted");
}

} // namespace

int main()
{
    try {
        test_full_rank();
        test_rank_deficient();
        test_no_measurements();
        test_no_feasible_directions();
        test_invalid_basis();

        std::cout
            << "PASS E5.5c.1 measurement information\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E5.5c.1: "
            << error.what() << '\n';

        return 1;
    }
}
