#include "tpsf/dense_full_space_sqp.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

using namespace tpsf;
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void test_unconstrained()
{
    const auto assemble = [](const Vector& u) {
        DenseFullSpaceSubproblem s;

        s.constraint_jacobian = Matrix(0, 2);
        s.constraint_residual = Vector(0);

        s.measurement_jacobian = Matrix::Identity(2, 2);
        s.measurement_residual = u;

        s.measurement_residual[0] -= 2.0;
        s.measurement_residual[1] += 3.0;

        return s;
    };

    Vector initial(2);
    initial << 10.0, 10.0;

    const auto result =
        solve_dense_full_space_sqp(assemble, initial);

    require(result.converged(),
            "Unconstrained problem did not converge");

    require((result.coordinates -
             (Vector(2) << 2.0, -3.0).finished()).norm()
                < 1e-10,
            "Incorrect unconstrained solution");
}

void test_nonlinear_equality()
{
    // Minimize 1/2 (x-2)^2 subject to x^2=1.
    // Starting near +1 selects the positive branch.
    const auto assemble = [](const Vector& u) {
        DenseFullSpaceSubproblem s;

        const double x = u[0];

        s.constraint_jacobian = Matrix(1, 1);
        s.constraint_jacobian(0, 0) = 2.0*x;

        s.constraint_residual = Vector(1);
        s.constraint_residual[0] = x*x-1.0;

        s.measurement_jacobian = Matrix(1, 1);
        s.measurement_jacobian(0, 0) = 1.0;

        s.measurement_residual = Vector(1);
        s.measurement_residual[0] = x-2.0;

        return s;
    };

    Vector initial(1);
    initial << 1.2;

    const auto result =
        solve_dense_full_space_sqp(assemble, initial);

    require(result.converged(),
            "Nonlinear equality did not converge");

    require(std::abs(result.coordinates[0]-1.0) < 1e-8,
            "Incorrect nonlinear equality solution");

    require(result.final_constraint_l2 < 1e-8,
            "Nonlinear equality remains infeasible");
}

void test_incompatible_linearization()
{
    const auto assemble = [](const Vector& u) {
        DenseFullSpaceSubproblem s;

        s.constraint_jacobian = Matrix(2, 1);
        s.constraint_jacobian << 1.0, 1.0;

        s.constraint_residual = Vector(2);
        s.constraint_residual << u[0]-1.0, u[0]-2.0;

        s.measurement_jacobian = Matrix(0, 1);
        s.measurement_residual = Vector(0);

        return s;
    };

    Vector initial = Vector::Zero(1);

    const auto result =
        solve_dense_full_space_sqp(assemble, initial);

    require(
        result.termination ==
            DenseSQPTermination::incompatible_linearization,
        "Incompatible linearization was accepted");
}

void test_unresolved_direction()
{
    const auto assemble = [](const Vector& u) {
        DenseFullSpaceSubproblem s;

        s.constraint_jacobian = Matrix(1, 3);
        s.constraint_jacobian << 1.0, 0.0, 0.0;

        s.constraint_residual = Vector(1);
        s.constraint_residual[0] = u[0]-2.0;

        s.measurement_jacobian = Matrix(1, 3);
        s.measurement_jacobian << 0.0, 1.0, 0.0;

        s.measurement_residual = Vector(1);
        s.measurement_residual[0] = u[1]-4.0;

        return s;
    };

    Vector initial = Vector::Zero(3);

    const auto result =
        solve_dense_full_space_sqp(assemble, initial);

    require(result.converged(),
            "Rank-deficient problem did not converge");

    require(std::abs(result.coordinates[0]-2.0) < 1e-10,
            "Incorrect constrained coordinate");

    require(std::abs(result.coordinates[1]-4.0) < 1e-10,
            "Incorrect observed coordinate");

    require(std::abs(result.coordinates[2]) < 1e-10,
            "Unresolved coordinate unexpectedly changed");

    require(!result.history.empty(),
            "Missing rank diagnostics");

    require(result.history.front().unresolved_directions == 1,
            "Unresolved direction was not reported");
}

void test_iteration_limit()
{
    const auto assemble = [](const Vector& u) {
        DenseFullSpaceSubproblem s;

        s.constraint_jacobian = Matrix(1, 1);
        s.constraint_jacobian(0, 0) = 2.0*u[0];

        s.constraint_residual = Vector(1);
        s.constraint_residual[0] = u[0]*u[0]-1.0;

        s.measurement_jacobian = Matrix(0, 1);
        s.measurement_residual = Vector(0);

        return s;
    };

    DenseSQPOptions options;
    options.maximum_iterations = 1;

    Vector initial(1);
    initial << 2.0;

    const auto result =
        solve_dense_full_space_sqp(
            assemble, initial, options);

    require(
        result.termination ==
            DenseSQPTermination::maximum_iterations,
        "Iteration limit was not enforced");
}


void test_zero_dimensional_feasible()
{
    const auto assemble = [](const Vector&) {
        DenseFullSpaceSubproblem s;

        s.constraint_jacobian = Matrix(1, 0);
        s.constraint_residual = Vector::Zero(1);

        s.measurement_jacobian = Matrix(0, 0);
        s.measurement_residual = Vector(0);

        return s;
    };

    const auto result =
        solve_dense_full_space_sqp(
            assemble, Vector(0));

    require(result.converged(),
            "Feasible zero-dimensional problem failed");
}

void test_zero_dimensional_infeasible()
{
    const auto assemble = [](const Vector&) {
        DenseFullSpaceSubproblem s;

        s.constraint_jacobian = Matrix(1, 0);
        s.constraint_residual = Vector::Ones(1);

        s.measurement_jacobian = Matrix(0, 0);
        s.measurement_residual = Vector(0);

        return s;
    };

    const auto result =
        solve_dense_full_space_sqp(
            assemble, Vector(0));

    require(
        result.termination ==
            DenseSQPTermination::incompatible_linearization,
        "Infeasible zero-dimensional problem accepted");
}

void test_nonlinear_rank_loss()
{
    // c(x)=x^2-1, starting at x=0.
    // The Jacobian vanishes while the constraint is
    // infeasible. This must not be reported as convergence.
    const auto assemble = [](const Vector& u) {
        DenseFullSpaceSubproblem s;

        s.constraint_jacobian = Matrix(1, 1);
        s.constraint_jacobian(0, 0) = 2.0*u[0];

        s.constraint_residual = Vector(1);
        s.constraint_residual[0] = u[0]*u[0]-1.0;

        s.measurement_jacobian = Matrix(0, 1);
        s.measurement_residual = Vector(0);

        return s;
    };

    const auto result =
        solve_dense_full_space_sqp(
            assemble, Vector::Zero(1));

    require(!result.converged(),
            "Rank loss produced false convergence");

    require(
        result.termination ==
            DenseSQPTermination::incompatible_linearization,
        "Rank-loss termination not reported");
}

void test_line_search_failure()
{
    // Deliberately inconsistent nonlinear function and
    // Jacobian. The linear model predicts improvement,
    // but every nonzero trial worsens the true objective.
    //
    // This verifies rejection behavior independently of
    // the exact derivative tests elsewhere.
    const auto assemble = [](const Vector& u) {
        DenseFullSpaceSubproblem s;

        s.constraint_jacobian = Matrix(0, 1);
        s.constraint_residual = Vector(0);

        s.measurement_jacobian = Matrix(1, 1);
        s.measurement_jacobian(0, 0) = 1.0;

        s.measurement_residual = Vector(1);
        s.measurement_residual[0] = 1.0 + u[0]*u[0];

        return s;
    };

    DenseSQPOptions options;
    options.maximum_backtracks = 8;

    Vector initial = Vector::Zero(1);

    const auto result =
        solve_dense_full_space_sqp(
            assemble, initial, options);

    require(
        result.termination ==
            DenseSQPTermination::line_search_failure,
        "Line-search failure was not reported");

    require(result.coordinates.norm() == 0.0,
            "Rejected trial modified accepted coordinates");
}

void test_nonlinear_backtracking()
{
    // Newton's method for x^2=1, initialized at x=0.2.
    // The full Newton step overshoots substantially.
    const auto assemble = [](const Vector& u) {
        DenseFullSpaceSubproblem s;

        const double x = u[0];

        s.constraint_jacobian = Matrix(1, 1);
        s.constraint_jacobian(0, 0) = 2.0*x;

        s.constraint_residual = Vector(1);
        s.constraint_residual[0] = x*x-1.0;

        s.measurement_jacobian = Matrix(0, 1);
        s.measurement_residual = Vector(0);

        return s;
    };

    Vector initial(1);
    initial << 0.2;

    const auto result =
        solve_dense_full_space_sqp(
            assemble, initial);

    require(result.converged(),
            "Backtracking problem did not converge");

    require(std::abs(result.coordinates[0]-1.0) < 1e-8,
            "Incorrect backtracking solution");

    bool backtracked = false;

    for (const auto& iteration : result.history) {
        if (iteration.step_length > 0.0 &&
            iteration.step_length < 1.0)
            backtracked = true;
    }

    require(backtracked,
            "Expected nonlinear backtracking did not occur");
}

} // namespace

int main()
{
    try {
        test_unconstrained();
        test_nonlinear_equality();
        test_incompatible_linearization();
        test_unresolved_direction();
        test_iteration_limit();
        test_zero_dimensional_feasible();
        test_zero_dimensional_infeasible();
        test_nonlinear_rank_loss();
        test_line_search_failure();
        test_nonlinear_backtracking();

        std::cout
            << "PASS E2.4b.11b.2 nonlinear SQP reference solver\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E2.4b.11b.2: "
            << error.what() << '\n';

        return 1;
    }
}
