#include "tpsf/dense_constrained_least_squares.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using namespace tpsf;
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(double a, double b, double tolerance,
          const std::string& message)
{
    require(std::isfinite(a) && std::isfinite(b),
            message + ": nonfinite");

    require(
        std::abs(a-b) <=
        tolerance*std::max({1.0, std::abs(a), std::abs(b)}),
        message);
}

void test_unique_solution()
{
    Matrix A(1, 2);
    A << 1.0, 1.0;

    Vector c(1);
    c << -3.0;

    Matrix B = Matrix::Identity(2, 2);
    Vector r = Vector::Zero(2);

    const auto result =
        solve_dense_constrained_least_squares(A, c, B, r);

    require(result.linearized_constraints_compatible,
            "Unique solution: incompatible");

    require(result.constraint_rank == 1,
            "Unique solution: constraint rank");

    require(result.reduced_measurement_rank == 1,
            "Unique solution: reduced rank");

    require(result.unresolved_directions == 0,
            "Unique solution: unresolved directions");

    near(result.step[0], 1.5, 1e-12, "Unique x");
    near(result.step[1], 1.5, 1e-12, "Unique y");
}

void test_redundant_constraints()
{
    Matrix A(2, 2);
    A << 1.0, 1.0,
         2.0, 2.0;

    Vector c(2);
    c << -3.0, -6.0;

    Matrix B = Matrix::Identity(2, 2);
    Vector r = Vector::Zero(2);

    const auto result =
        solve_dense_constrained_least_squares(A, c, B, r);

    require(result.linearized_constraints_compatible,
            "Redundant constraints: incompatible");

    require(result.constraint_rank == 1,
            "Redundant constraints: incorrect rank");

    near(result.step[0], 1.5, 1e-12, "Redundant x");
    near(result.step[1], 1.5, 1e-12, "Redundant y");
}

void test_inconsistent_constraints()
{
    Matrix A(2, 1);
    A << 1.0,
         1.0;

    Vector c(2);
    c << -1.0, -2.0;

    Matrix B(0, 1);
    Vector r(0);

    const auto result =
        solve_dense_constrained_least_squares(A, c, B, r);

    require(!result.linearized_constraints_compatible,
            "Inconsistent constraints accepted");

    require(result.linearized_constraint_residual > 0.1,
            "Inconsistent residual not reported");
}

void test_unobservable_direction()
{
    Matrix A(1, 3);
    A << 1.0, 0.0, 0.0;

    Vector c(1);
    c << -2.0;

    Matrix B(1, 3);
    B << 0.0, 1.0, 0.0;

    Vector r(1);
    r << -4.0;

    const auto result =
        solve_dense_constrained_least_squares(A, c, B, r);

    require(result.linearized_constraints_compatible,
            "Unobservable case incompatible");

    require(result.constraint_rank == 1,
            "Unobservable constraint rank");

    require(result.reduced_measurement_rank == 1,
            "Unobservable measurement rank");

    require(result.unresolved_directions == 1,
            "Unobservable direction not reported");

    near(result.step[0], 2.0, 1e-12, "Unobservable x");
    near(result.step[1], 4.0, 1e-12, "Unobservable y");
    near(result.step[2], 0.0, 1e-12, "Unobservable z");
}

void test_no_measurements()
{
    Matrix A(1, 2);
    A << 1.0, 1.0;

    Vector c(1);
    c << -4.0;

    Matrix B(0, 2);
    Vector r(0);

    const auto result =
        solve_dense_constrained_least_squares(A, c, B, r);

    require(result.linearized_constraints_compatible,
            "No-measurement case incompatible");

    require(result.reduced_measurement_rank == 0,
            "No-measurement rank");

    require(result.unresolved_directions == 1,
            "No-measurement nullity");

    near(result.step[0], 2.0, 1e-12, "No-measurement x");
    near(result.step[1], 2.0, 1e-12, "No-measurement y");
}

void test_singular_measurements()
{
    Matrix A(0, 2);
    Vector c(0);

    Matrix B(2, 2);
    B << 1.0, 0.0,
         2.0, 0.0;

    Vector r(2);
    r << -3.0, -6.0;

    const auto result =
        solve_dense_constrained_least_squares(A, c, B, r);

    require(result.linearized_constraints_compatible,
            "Singular measurement case incompatible");

    require(result.constraint_rank == 0,
            "Singular measurement constraint rank");

    require(result.reduced_measurement_rank == 1,
            "Singular measurement reduced rank");

    require(result.unresolved_directions == 1,
            "Singular measurement nullity");

    near(result.step[0], 3.0, 1e-12, "Singular x");
    near(result.step[1], 0.0, 1e-12, "Singular y");
}

void test_invalid_dimensions()
{
    Matrix A(1, 2);
    A.setZero();

    Vector c(2);
    c.setZero();

    Matrix B(0, 2);
    Vector r(0);

    bool caught = false;

    try {
        solve_dense_constrained_least_squares(A, c, B, r);
    }
    catch (const std::invalid_argument&) {
        caught = true;
    }

    require(caught, "Invalid dimensions accepted");
}

} // namespace

int main()
{
    try {
        test_unique_solution();
        test_redundant_constraints();
        test_inconsistent_constraints();
        test_unobservable_direction();
        test_no_measurements();
        test_singular_measurements();
        test_invalid_dimensions();

        std::cout
            << "PASS E2.4b.11a dense constrained least squares\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E2.4b.11a: "
            << error.what() << '\n';

        return 1;
    }
}
