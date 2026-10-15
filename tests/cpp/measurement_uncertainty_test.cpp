#include "tpsf/measurement_uncertainty.hpp"

#include <Eigen/Dense>

#include <algorithm>
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

void near(double actual,
          double expected,
          double tolerance,
          const std::string& message)
{
    require(
        std::abs(actual-expected) <= tolerance,
        message);
}

void test_full_rank()
{
    // Two independent measurements:
    //
    // S = diag(2,4)
    // Cov(z) = diag(1/4,1/16)
    //
    // q = 3*z1 + 5*z2
    //
    // Var(q) = 9/4 + 25/16 = 3.8125

    Eigen::Matrix2d S = Eigen::Matrix2d::Zero();
    S(0,0) = 2.0;
    S(1,1) = 4.0;

    Eigen::Vector2d g;
    g << 3.0, 5.0;

    const auto result =
        tpsf::propagate_measurement_uncertainty(S, g);

    require(result.status ==
                tpsf::UncertaintyStatus::finite,
            "Full-rank case should be finite");

    require(result.information_rank == 2,
            "Incorrect full-rank information rank");

    near(result.variance, 3.8125, 1.0e-12,
         "Incorrect analytic variance");

    near(result.standard_deviation,
         std::sqrt(3.8125), 1.0e-12,
         "Incorrect analytic standard deviation");
}

void test_unresolved_quantity()
{
    Eigen::Matrix<double,1,2> S;
    S << 2.0, 0.0;

    Eigen::Vector2d g;
    g << 3.0, 5.0;

    const auto result =
        tpsf::propagate_measurement_uncertainty(S, g);

    require(result.status ==
                tpsf::UncertaintyStatus::unbounded,
            "Unresolved quantity must be unbounded");

    require(result.unresolved_directions == 1,
            "Expected one unresolved direction");

    near(result.unresolved_sensitivity_norm,
         5.0, 1.0e-12,
         "Incorrect unresolved sensitivity");

    require(std::isinf(result.variance),
            "Unbounded variance must be infinite");
}

void test_rank_deficient_but_quantity_resolved()
{
    Eigen::Matrix<double,1,2> S;
    S << 2.0, 0.0;

    Eigen::Vector2d g;
    g << 3.0, 0.0;

    const auto result =
        tpsf::propagate_measurement_uncertainty(S, g);

    require(result.status ==
                tpsf::UncertaintyStatus::finite,
            "Quantity should be identifiable");

    require(result.unresolved_directions == 1,
            "State should remain incompletely observable");

    near(result.variance, 2.25, 1.0e-12,
         "Incorrect rank-deficient analytic variance");
}

void test_no_measurements()
{
    Eigen::MatrixXd S(0,2);

    Eigen::Vector2d g;
    g << 1.0, 0.0;

    const auto result =
        tpsf::propagate_measurement_uncertainty(S, g);

    require(result.status ==
                tpsf::UncertaintyStatus::unbounded,
            "No measurements must leave quantity unbounded");

    g.setZero();

    const auto fixed =
        tpsf::propagate_measurement_uncertainty(S, g);

    require(fixed.status ==
                tpsf::UncertaintyStatus::finite,
            "Locally insensitive quantity should be finite");

    near(fixed.variance, 0.0, 0.0,
         "Locally fixed quantity must have zero variance");
}

void test_no_feasible_directions()
{
    Eigen::MatrixXd S(3,0);
    Eigen::VectorXd g(0);

    const auto result =
        tpsf::propagate_measurement_uncertainty(S, g);

    require(result.status ==
                tpsf::UncertaintyStatus::finite,
            "Constrained quantity should be locally fixed");

    near(result.variance, 0.0, 0.0,
         "Zero-dimensional variance must be zero");
}

void test_invalid_input()
{
    Eigen::MatrixXd S(1,2);
    S.setZero();

    Eigen::VectorXd g(3);
    g.setZero();

    bool rejected = false;

    try {
        (void)tpsf::propagate_measurement_uncertainty(S, g);
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    require(rejected, "Dimension mismatch was not rejected");
}

} // namespace

int main()
{
    try {
        test_full_rank();
        test_unresolved_quantity();
        test_rank_deficient_but_quantity_resolved();
        test_no_measurements();
        test_no_feasible_directions();
        test_invalid_input();

        std::cout
            << "PASS E5.3 measurement uncertainty propagation\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E5.3: " << error.what() << '\n';

        return 1;
    }
}
