#include "tpsf/matrix_free_kkt.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace tpsf {
namespace {

void validate_finite(const std::vector<double>& values)
{
    for (double value : values)
        if (!std::isfinite(value))
            throw std::invalid_argument(
                "Nonfinite KKT input or result");
}

} // namespace

MatrixFreeKKT::MatrixFreeKKT(
    const FullSpaceProblem& problem,
    std::vector<double> coordinates)
    : problem_(problem),
      coordinates_(std::move(coordinates))
{
    (void)problem_.trajectory(coordinates_);
}

std::size_t MatrixFreeKKT::primal_size() const noexcept
{
    return problem_.free_variable_count();
}

std::size_t MatrixFreeKKT::dual_size() const noexcept
{
    return problem_.constraint_count();
}

std::vector<double> MatrixFreeKKT::apply(
    const std::vector<double>& vector) const
{
    const auto n = primal_size();
    const auto m = dual_size();

    if (vector.size() != n+m)
        throw std::invalid_argument(
            "Incorrect KKT vector dimension");

    validate_finite(vector);

    const std::vector<double> primal(
        vector.begin(), vector.begin()+n);

    const std::vector<double> dual(
        vector.begin()+n, vector.end());

    const auto Bv =
        problem_.measurement_jvp(coordinates_, primal);

    auto top =
        problem_.measurement_vjp(coordinates_, Bv);

    const auto Atw =
        problem_.constraint_vjp(coordinates_, dual);

    const auto bottom =
        problem_.constraint_jvp(coordinates_, primal);

    if (top.size() != n ||
        Atw.size() != n ||
        bottom.size() != m)
        throw std::runtime_error(
            "Inconsistent KKT operator dimensions");

    for (std::size_t i = 0; i < n; ++i)
        top[i] += Atw[i];

    top.insert(top.end(), bottom.begin(), bottom.end());

    validate_finite(top);

    return top;
}

std::vector<double> MatrixFreeKKT::right_hand_side() const
{
    const auto residual =
        problem_.measurement_residual(coordinates_);

    const auto evaluation =
        problem_.evaluate(coordinates_);

    auto result =
        problem_.measurement_vjp(coordinates_, residual);

    if (result.size() != primal_size() ||
        evaluation.constraints.size() != dual_size())
        throw std::runtime_error(
            "Inconsistent KKT right-hand-side dimensions");

    for (double& value : result)
        value = -value;

    for (double value : evaluation.constraints)
        result.push_back(-value);

    validate_finite(result);

    return result;
}

} // namespace tpsf
