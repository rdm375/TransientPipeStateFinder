#include "tpsf/engineering_sensitivity.hpp"
#include "pipe_sim/integrate.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace tpsf {
namespace {

std::vector<double> to_std(const Eigen::VectorXd& u)
{
    return std::vector<double>(
        u.data(), u.data() + u.size());
}

std::vector<double> pressures(const std::vector<double>& state,
                              int cells)
{
    std::vector<double> p(
        static_cast<std::size_t>(cells + 1));

    p[0] = state[0];

    for (int i = 0; i < cells; ++i)
        p[static_cast<std::size_t>(i + 1)] =
            state[static_cast<std::size_t>(2*i + 3)];

    return p;
}

struct QuantityValue {
    double value = 0.0;
    std::size_t active_count = 0;
};

QuantityValue evaluate(
    const FullSpaceProblem& problem,
    const pipe_sim::TransientParameters& parameters,
    const Eigen::VectorXd& coordinates,
    EngineeringQuantity quantity,
    double tie_tolerance)
{
    if (coordinates.size() !=
            static_cast<Eigen::Index>(
                problem.free_variable_count()) ||
        !coordinates.allFinite())
        throw std::invalid_argument(
            "Invalid engineering sensitivity coordinates");

    const auto trajectory =
        problem.trajectory(to_std(coordinates));

    if (trajectory.states.empty())
        throw std::runtime_error("Empty trajectory");

    if (quantity == EngineeringQuantity::minimum_pressure) {
        double minimum = std::numeric_limits<double>::infinity();

        for (const auto& state : trajectory.states) {
            for (double p : pressures(state, trajectory.n))
                minimum = std::min(minimum, p);
        }

        std::size_t active = 0;

        for (const auto& state : trajectory.states) {
            for (double p : pressures(state, trajectory.n)) {
                if (p - minimum <= tie_tolerance)
                    ++active;
            }
        }

        if (!std::isfinite(minimum))
            throw std::runtime_error(
                "Nonfinite minimum pressure");

        return {minimum, active};
    }

    const auto& state =
        quantity == EngineeringQuantity::initial_linepack
            ? trajectory.states.front()
            : trajectory.states.back();

    const double mass = pipe_sim::linepack_cz(
        trajectory.n,
        pressures(state, trajectory.n),
        parameters);

    if (!std::isfinite(mass))
        throw std::runtime_error("Nonfinite linepack");

    return {mass, 0};
}

} // namespace

double evaluate_engineering_quantity(
    const FullSpaceProblem& problem,
    const pipe_sim::TransientParameters& parameters,
    const Eigen::VectorXd& coordinates,
    EngineeringQuantity quantity)
{
    return evaluate(problem, parameters, coordinates, quantity, 0.0).value;
}

QuantitySensitivity assess_engineering_sensitivity(
    const FullSpaceProblem& problem,
    const pipe_sim::TransientParameters& parameters,
    const Eigen::VectorXd& coordinates,
    const StatefindingQuality& quality,
    EngineeringQuantity quantity,
    const EngineeringSensitivityOptions& options)
{
    if (!std::isfinite(options.difference_step) ||
        options.difference_step <= 0.0 ||
        !std::isfinite(options.minimum_pressure_tie_tolerance_pa) ||
        options.minimum_pressure_tie_tolerance_pa < 0.0)
        throw std::invalid_argument(
            "Invalid sensitivity options");

    const Eigen::Index n = coordinates.size();

    if (quality.feasible_basis.rows() != n ||
        quality.unresolved_basis.rows() != n)
        throw std::invalid_argument(
            "Observability basis dimension mismatch");

    QuantitySensitivity result;
    result.quantity = quantity;

    const auto nominal = evaluate(
        problem,
        parameters,
        coordinates,
        quantity,
        options.minimum_pressure_tie_tolerance_pa);

    result.value = nominal.value;
    result.active_minimum_count = nominal.active_count;

    if (quantity == EngineeringQuantity::minimum_pressure &&
        nominal.active_count != 1) {
        result.differentiable = false;
        return result;
    }

    result.gradient = Eigen::VectorXd::Zero(n);

    for (Eigen::Index j = 0; j < n; ++j) {
        const double h =
            options.difference_step *
            std::max(1.0, std::abs(coordinates[j]));

        Eigen::VectorXd plus = coordinates;
        Eigen::VectorXd minus = coordinates;

        plus[j] += h;
        minus[j] -= h;

        const double qp =
            evaluate(problem, parameters, plus, quantity, 0.0).value;

        const double qm =
            evaluate(problem, parameters, minus, quantity, 0.0).value;

        result.gradient[j] = (qp - qm) / (2.0*h);
    }

    result.feasible_sensitivity =
        quality.feasible_basis.transpose() * result.gradient;

    result.unresolved_sensitivity =
        quality.unresolved_basis.transpose() * result.gradient;

    result.feasible_sensitivity_norm =
        result.feasible_sensitivity.norm();

    result.unresolved_sensitivity_norm =
        result.unresolved_sensitivity.norm();

    if (!result.gradient.allFinite() ||
        !std::isfinite(result.feasible_sensitivity_norm) ||
        !std::isfinite(result.unresolved_sensitivity_norm))
        throw std::runtime_error(
            "Nonfinite engineering sensitivity");

    return result;
}

} // namespace tpsf
