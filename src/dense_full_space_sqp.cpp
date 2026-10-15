#include "tpsf/dense_full_space_sqp.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace tpsf {
namespace {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

void validate_options(const DenseSQPOptions& o)
{
    if (o.maximum_iterations == 0 ||
        o.maximum_backtracks == 0 ||
        !std::isfinite(o.feasibility_tolerance) ||
        o.feasibility_tolerance < 0.0 ||
        !std::isfinite(o.stationarity_tolerance) ||
        o.stationarity_tolerance < 0.0 ||
        !std::isfinite(o.initial_penalty) ||
        o.initial_penalty <= 0.0 ||
        !std::isfinite(o.penalty_margin) ||
        o.penalty_margin <= 0.0 ||
        !std::isfinite(o.backtracking_factor) ||
        o.backtracking_factor <= 0.0 ||
        o.backtracking_factor >= 1.0 ||
        !std::isfinite(o.armijo_fraction) ||
        o.armijo_fraction <= 0.0 ||
        o.armijo_fraction >= 1.0)
        throw std::invalid_argument("Invalid SQP options");
}

void validate(const DenseFullSpaceSubproblem& s,
              Eigen::Index n)
{
    if (s.constraint_jacobian.cols() != n ||
        s.measurement_jacobian.cols() != n ||
        s.constraint_jacobian.rows() !=
            s.constraint_residual.size() ||
        s.measurement_jacobian.rows() !=
            s.measurement_residual.size() ||
        !s.constraint_jacobian.allFinite() ||
        !s.measurement_jacobian.allFinite() ||
        !s.constraint_residual.allFinite() ||
        !s.measurement_residual.allFinite())
        throw std::runtime_error(
            "Invalid nonlinear SQP assembly");
}

double objective(const DenseFullSpaceSubproblem& s)
{
    return 0.5*s.measurement_residual.squaredNorm();
}

double constraint_l1(const DenseFullSpaceSubproblem& s)
{
    return s.constraint_residual.lpNorm<1>();
}

double merit(const DenseFullSpaceSubproblem& s,
             double penalty)
{
    return objective(s)+penalty*constraint_l1(s);
}

// Compute the component of the objective gradient
// orthogonal to the numerical range of A^T.
double projected_gradient_norm(
    const DenseFullSpaceSubproblem& s,
    double rank_tolerance)
{
    const Vector gradient =
        s.measurement_jacobian.transpose() *
        s.measurement_residual;

    if (gradient.size() == 0 ||
        s.constraint_jacobian.rows() == 0)
        return gradient.norm();

    const Matrix transpose =
        s.constraint_jacobian.transpose();

    Eigen::JacobiSVD<Matrix> svd(
        transpose,
        Eigen::ComputeThinU | Eigen::ComputeThinV);

    const auto singular = svd.singularValues();

    Eigen::Index rank = 0;

    if (singular.size() > 0 && singular[0] > 0.0) {
        for (Eigen::Index i = 0; i < singular.size(); ++i)
            if (singular[i] >
                rank_tolerance*singular[0])
                ++rank;
    }

    Vector projected = gradient;

    if (rank > 0) {
        const auto basis = svd.matrixU().leftCols(rank);
        projected -= basis*(basis.transpose()*gradient);
    }

    return projected.norm();
}


double constraint_l1_directional_derivative(
    const Vector& residual,
    const Vector& directional_change)
{
    if (residual.size() != directional_change.size())
        throw std::invalid_argument(
            "Incompatible directional derivative dimensions");

    double derivative = 0.0;

    for (Eigen::Index i = 0; i < residual.size(); ++i) {
        if (residual[i] > 0.0)
            derivative += directional_change[i];
        else if (residual[i] < 0.0)
            derivative -= directional_change[i];
        else
            derivative += std::abs(directional_change[i]);
    }

    return derivative;
}

DenseSQPIteration diagnostics(
    const DenseFullSpaceSubproblem& s,
    std::size_t iteration,
    double projected_gradient,
    double penalty)
{
    DenseSQPIteration d;
    d.iteration = iteration;
    d.objective = objective(s);
    d.constraint_l1 = constraint_l1(s);
    d.constraint_l2 = s.constraint_residual.norm();
    d.projected_gradient = projected_gradient;
    d.penalty = penalty;
    return d;
}

} // namespace

DenseSQPResult solve_dense_full_space_sqp(
    const DenseSQPAssembly& assemble,
    const Vector& initial,
    const DenseSQPOptions& options)
{
    validate_options(options);

    if (!assemble || !initial.allFinite())
        throw std::invalid_argument(
            "Invalid SQP callback or initial coordinates");

    DenseSQPResult result;
    result.coordinates = initial;

    double penalty = options.initial_penalty;

    auto current = assemble(result.coordinates);
    validate(current, initial.size());

    for (std::size_t k = 0;
         k <= options.maximum_iterations; ++k) {

        const double stationarity =
            projected_gradient_norm(
                current,
                options.linear_solver.relative_rank_tolerance);

        auto record = diagnostics(
            current, k, stationarity, penalty);

        if (!std::isfinite(record.objective) ||
            !std::isfinite(record.constraint_l1) ||
            !std::isfinite(record.constraint_l2) ||
            !std::isfinite(stationarity)) {
            result.termination =
                DenseSQPTermination::numerical_failure;
            result.history.push_back(record);
            break;
        }

        const Vector gradient =
            current.measurement_jacobian.transpose() *
            current.measurement_residual;

        const bool feasible =
            record.constraint_l2 <=
            options.feasibility_tolerance;

        const bool stationary =
            stationarity <=
            options.stationarity_tolerance *
            std::max(1.0, gradient.norm());

        if (feasible && stationary) {
            result.termination =
                DenseSQPTermination::converged;
            result.history.push_back(record);
            break;
        }

        if (k == options.maximum_iterations) {
            result.termination =
                DenseSQPTermination::maximum_iterations;
            result.history.push_back(record);
            break;
        }

        const auto linear =
            solve_dense_constrained_least_squares(
                current.constraint_jacobian,
                current.constraint_residual,
                current.measurement_jacobian,
                current.measurement_residual,
                options.linear_solver);

        record.constraint_rank = linear.constraint_rank;
        record.reduced_measurement_rank =
            linear.reduced_measurement_rank;
        record.unresolved_directions =
            linear.unresolved_directions;

        if (!linear.linearized_constraints_compatible) {
            result.termination =
                DenseSQPTermination::incompatible_linearization;
            result.history.push_back(record);
            break;
        }

        const Vector& step = linear.step;
        record.step_norm = step.norm();

        if (!step.allFinite() ||
            !std::isfinite(record.step_norm)) {
            result.termination =
                DenseSQPTermination::numerical_failure;
            result.history.push_back(record);
            break;
        }

        const Vector Ad = current.constraint_jacobian*step;
        const Vector Bd = current.measurement_jacobian*step;

        const double objective_derivative =
            gradient.dot(step);

        const double constraint_derivative =
            constraint_l1_directional_derivative(
                current.constraint_residual, Ad);

        // The directional derivative of the exact L1 merit
        // function is:
        //
        // D Phi = grad(J)^T d + penalty * D ||c||_1.
        //
        // Increase the penalty only when the constraint
        // contribution is strictly decreasing.
        if (constraint_derivative < 0.0) {
            const double margin =
                options.penalty_margin *
                std::max(1.0, std::abs(objective_derivative));

            const double required =
                (objective_derivative + margin) /
                (-constraint_derivative);

            if (required > penalty)
                penalty = required;
        }

        const double merit_derivative =
            objective_derivative +
            penalty*constraint_derivative;

        if (!std::isfinite(penalty) ||
            !std::isfinite(merit_derivative)) {
            result.termination =
                DenseSQPTermination::numerical_failure;
            result.history.push_back(record);
            break;
        }

        // A zero or positive directional derivative does not
        // establish a globalization direction.
        if (!(merit_derivative < 0.0)) {
            result.termination =
                DenseSQPTermination::non_descent_direction;
            result.history.push_back(record);
            break;
        }

        record.penalty = penalty;

        const double current_merit = merit(current, penalty);

        bool accepted = false;
        double alpha = 1.0;

        for (std::size_t trial = 0;
             trial < options.maximum_backtracks; ++trial) {

            const double armijo_bound =
                current_merit +
                options.armijo_fraction *
                alpha*merit_derivative;

            const Vector candidate =
                result.coordinates+alpha*step;

            if (!candidate.allFinite()) {
                alpha *= options.backtracking_factor;
                continue;
            }

            // A failed physical trial is rejected, not
            // silently substituted with an artificial value.
            try {
                auto trial_problem = assemble(candidate);
                validate(trial_problem, initial.size());

                const double actual_merit =
                    merit(trial_problem, penalty);

                if (std::isfinite(actual_merit) &&
                    actual_merit <= armijo_bound) {

                    result.coordinates = candidate;
                    current = std::move(trial_problem);

                    record.step_length = alpha;
                    accepted = true;
                    break;
                }
            }
            catch (const std::exception&) {
                // Reject this trial and backtrack.
            }

            alpha *= options.backtracking_factor;
        }

        result.history.push_back(record);

        if (!accepted) {
            result.termination =
                DenseSQPTermination::line_search_failure;
            break;
        }
    }

    result.final_objective = objective(current);
    result.final_constraint_l2 =
        current.constraint_residual.norm();

    result.final_projected_gradient =
        projected_gradient_norm(
            current,
            options.linear_solver.relative_rank_tolerance);

    return result;
}

DenseSQPResult solve_dense_full_space_sqp(
    const FullSpaceProblem& problem,
    const Vector& initial,
    const DenseSQPOptions& options)
{
    return solve_dense_full_space_sqp(
        [&problem](const Vector& u) {
            return assemble_dense_full_space_subproblem(
                problem, u);
        },
        initial,
        options);
}

} // namespace tpsf
