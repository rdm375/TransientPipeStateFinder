#include "tpsf/minres.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

tpsf::SymmetricOperator operator_from_matrix(
    const Matrix& matrix)
{
    return [matrix](const tpsf::KrylovVector& x) {
        if (static_cast<Eigen::Index>(x.size()) !=
            matrix.cols()) {
            throw std::invalid_argument("Incorrect input size");
        }

        Eigen::Map<const Vector> input(
            x.data(), static_cast<Eigen::Index>(x.size()));

        const Vector output = matrix*input;

        return tpsf::KrylovVector(
            output.data(), output.data()+output.size());
    };
}

void verify(
    const std::string& name,
    const Matrix& matrix,
    const Vector& rhs,
    bool expect_convergence)
{
    const auto apply = operator_from_matrix(matrix);

    const tpsf::KrylovVector b(
        rhs.data(), rhs.data()+rhs.size());

    tpsf::MinresOptions options;
    options.max_iterations = 100;
    options.relative_tolerance = 1e-10;
    options.absolute_tolerance = 1e-12;

    const auto result =
        tpsf::solve_minres(apply, b, options);

    const Eigen::Map<const Vector> x(
        result.solution.data(),
        static_cast<Eigen::Index>(result.solution.size()));

    const double residual = (matrix*x-rhs).norm();

    std::cout
        << name
        << " iterations=" << result.iterations
        << " residual=" << residual
        << " status=" << static_cast<int>(result.status)
        << '\n';

    if (expect_convergence) {
        require(
            result.status == tpsf::MinresStatus::converged,
            name + ": failed to converge");

        require(
            residual <= 1e-9*std::max(1.0, rhs.norm()),
            name + ": excessive residual");

        const Vector reference =
            matrix.completeOrthogonalDecomposition().solve(rhs);

        require(
            (x-reference).norm() <=
                1e-8*std::max(1.0, reference.norm()),
            name + ": dense reference mismatch");
    }
    else {
        require(
            result.status != tpsf::MinresStatus::converged,
            name + ": incorrectly reported convergence");
    }
}

} // namespace

int main()
{
    try {
        {
            Matrix A(3,3);
            A << 4,1,0,
                 1,3,1,
                 0,1,2;

            Vector b(3);
            b << 1,2,3;

            verify("SPD", A, b, true);
        }

        {
            Matrix A(3,3);
            A << 0,2,0,
                 2,0,1,
                 0,1,-3;

            Vector b(3);
            b << 1,-2,3;

            verify("indefinite", A, b, true);
        }

        {
            Matrix A = Matrix::Zero(3,3);
            A(0,0) = 2;
            A(1,1) = -1;

            Vector b(3);
            b << 4,-3,0;

            verify("singular-compatible", A, b, true);
        }

        {
            Matrix A = Matrix::Zero(3,3);
            A(0,0) = 2;
            A(1,1) = -1;

            Vector b(3);
            b << 4,-3,1;

            verify("singular-incompatible", A, b, false);
        }

        {
            const auto identity =
                [](const tpsf::KrylovVector& x) {
                    return x;
                };

            const auto result =
                tpsf::solve_minres(identity, {0,0,0});

            require(
                result.status == tpsf::MinresStatus::converged,
                "Zero RHS failed");
        }

        std::cout << "PASS E4.3b.1 MINRES tests\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E4.3b.1: "
            << error.what() << '\n';

        return 1;
    }
}
