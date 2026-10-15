#include "tpsf/dense_full_space_subproblem.hpp"

#include <stdexcept>
#include <vector>

namespace tpsf {
namespace {

Eigen::VectorXd to_eigen(const std::vector<double>& values)
{
    Eigen::VectorXd result(values.size());

    for (std::size_t i = 0; i < values.size(); ++i)
        result[static_cast<Eigen::Index>(i)] = values[i];

    return result;
}

std::vector<double> to_std(const Eigen::VectorXd& values)
{
    return std::vector<double>(
        values.data(), values.data()+values.size());
}

} // namespace

DenseFullSpaceSubproblem assemble_dense_full_space_subproblem(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& u)
{
    const auto n = static_cast<Eigen::Index>(
        problem.free_variable_count());

    const auto m = static_cast<Eigen::Index>(
        problem.constraint_count());

    if (u.size() != n || !u.allFinite())
        throw std::invalid_argument(
            "Invalid dense assembly coordinates");

    const auto coordinates = to_std(u);

    const auto evaluation = problem.evaluate(coordinates);

    DenseFullSpaceSubproblem result;

    result.constraint_residual =
        to_eigen(evaluation.constraints);

    result.measurement_residual =
        to_eigen(problem.measurement_residual(coordinates));

    const auto q = result.measurement_residual.size();

    result.constraint_jacobian.resize(m, n);
    result.measurement_jacobian.resize(q, n);

    std::vector<double> direction(
        static_cast<std::size_t>(n), 0.0);

    for (Eigen::Index j = 0; j < n; ++j) {
        direction[static_cast<std::size_t>(j)] = 1.0;

        const auto constraint_column =
            problem.constraint_jvp(coordinates, direction);

        const auto measurement_column =
            problem.measurement_jvp(coordinates, direction);

        if (static_cast<Eigen::Index>(
                constraint_column.size()) != m ||
            static_cast<Eigen::Index>(
                measurement_column.size()) != q)
            throw std::runtime_error(
                "Inconsistent assembled Jacobian dimensions");

        result.constraint_jacobian.col(j) =
            to_eigen(constraint_column);

        result.measurement_jacobian.col(j) =
            to_eigen(measurement_column);

        direction[static_cast<std::size_t>(j)] = 0.0;
    }

    if (!result.constraint_residual.allFinite() ||
        !result.measurement_residual.allFinite() ||
        !result.constraint_jacobian.allFinite() ||
        !result.measurement_jacobian.allFinite())
        throw std::runtime_error(
            "Nonfinite dense subproblem assembly");

    return result;
}

} // namespace tpsf
