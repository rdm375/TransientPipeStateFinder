#include "tpsf/objective.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace tpsf {

namespace {

double channel_cost(const std::vector<double>& residual)
{
    double cost = 0.0;

    for (double value : residual) {
        if (!std::isfinite(value))
            throw std::runtime_error("Nonfinite objective residual");

        cost += 0.5 * value * value;
    }

    if (!std::isfinite(cost))
        throw std::runtime_error("Nonfinite channel cost");

    return cost;
}

void validate_problem(const ObjectiveProblem& problem)
{
    const auto& cfg = problem.configuration;

    if (cfg.cells < 2 || cfg.steps < 0)
        throw std::invalid_argument("Invalid objective dimensions");

    validate_boundary_observations(problem.observations);

    const std::size_t count =
        static_cast<std::size_t>(cfg.steps) + 1;

    if (problem.observations.time.size() != count)
        throw std::invalid_argument(
            "Observation count does not match integration steps");

    if (!std::isfinite(cfg.parameters.dt) ||
        cfg.parameters.dt <= 0.0)
        throw std::invalid_argument("Invalid integration timestep");

    // E2.4a requires observations on the integration grid.
    // Irregular sampling will be supported separately.
    for (std::size_t k = 0; k < count; ++k) {
        const double expected =
            static_cast<double>(k) * cfg.parameters.dt;

        const double tolerance =
            1.0e-10 * std::max(1.0, std::abs(expected));

        if (std::abs(problem.observations.time[k] - expected)
            > tolerance)
            throw std::invalid_argument(
                "Observation timestamps differ from integration grid");
    }

    const auto& initial = problem.baseline_initial;

    if (initial.pressure.size() !=
            static_cast<std::size_t>(cfg.cells + 1) ||
        initial.flow.size() !=
            static_cast<std::size_t>(cfg.cells))
        throw std::invalid_argument(
            "Invalid baseline initial-state dimensions");

    if (!std::isfinite(initial.inlet_flow))
        throw std::invalid_argument(
            "Nonfinite baseline inlet flow");

    for (double value : initial.pressure) {
        if (!std::isfinite(value) || value <= 0.0)
            throw std::invalid_argument(
                "Invalid baseline initial pressure");
    }

    for (double value : initial.flow) {
        if (!std::isfinite(value))
            throw std::invalid_argument(
                "Nonfinite baseline initial flow");
    }

    if (!std::isfinite(problem.state_scales.pressure) ||
        !std::isfinite(problem.state_scales.flow) ||
        problem.state_scales.pressure <= 0.0 ||
        problem.state_scales.flow <= 0.0)
        throw std::invalid_argument("Invalid initial-state scales");
}

} // namespace

ObjectiveLayout objective_layout(
    const Configuration& cfg)
{
    if (cfg.cells < 2 || cfg.steps < 0)
        throw std::invalid_argument("Invalid objective dimensions");

    ObjectiveLayout layout;

    layout.initial_size =
        2 * static_cast<std::size_t>(cfg.cells) + 1;

    layout.boundary_size =
        static_cast<std::size_t>(cfg.steps) + 1;

    layout.inlet_pressure_offset = layout.initial_size;

    layout.outlet_flow_offset =
        layout.inlet_pressure_offset + layout.boundary_size;

    layout.total_size =
        layout.outlet_flow_offset + layout.boundary_size;

    return layout;
}

std::vector<double> zero_objective_parameters(
    const ObjectiveProblem& problem)
{
    validate_problem(problem);

    return std::vector<double>(
        objective_layout(problem.configuration).total_size,
        0.0);
}

ObjectiveEvaluation evaluate_objective(
    const ObjectiveProblem& problem,
    const std::vector<double>& parameters)
{
    validate_problem(problem);

    const auto layout =
        objective_layout(problem.configuration);

    if (parameters.size() != layout.total_size)
        throw std::invalid_argument(
            "Incorrect optimization-vector dimension");

    for (double value : parameters) {
        if (!std::isfinite(value))
            throw std::invalid_argument(
                "Nonfinite optimization parameter");
    }

    // Recover the scaled initial-state component.
    const std::vector<double> scaled_initial(
        parameters.begin(),
        parameters.begin() + layout.initial_size);

    InitialState initial = unpack_scaled_state(
        scaled_initial,
        problem.baseline_initial,
        problem.state_scales);

    BoundarySchedule schedule;

    schedule.inlet_pressure.resize(layout.boundary_size);
    schedule.outlet_flow.resize(layout.boundary_size);

    const auto& measured = problem.observations.measured;
    const auto& scales = problem.observations.scales;

    for (std::size_t k = 0; k < layout.boundary_size; ++k) {
        const double alpha =
            parameters[layout.inlet_pressure_offset + k];

        const double beta =
            parameters[layout.outlet_flow_offset + k];

        schedule.inlet_pressure[k] =
            measured.inlet_pressure[k] +
            scales.inlet_pressure * alpha;

        schedule.outlet_flow[k] =
            measured.outlet_flow[k] +
            scales.outlet_flow * beta;
    }

    // The prescribed initial inlet pressure is not an
    // independent initial-state optimization variable.
    initial.pressure[0] = schedule.inlet_pressure[0];

    ObjectiveEvaluation evaluation;

    evaluation.forward = run_forward_scheduled(
        problem.configuration,
        initial,
        schedule);

    if (evaluation.forward.info != 0 ||
        evaluation.forward.completed_steps !=
            problem.configuration.steps)
        throw std::runtime_error(
            "Forward integration failed during objective evaluation");

    evaluation.predicted.inlet_pressure =
        evaluation.forward.inlet_pressure;

    evaluation.predicted.inlet_flow =
        evaluation.forward.inlet_flow;

    evaluation.predicted.outlet_pressure =
        evaluation.forward.outlet_pressure;

    evaluation.predicted.outlet_flow =
        evaluation.forward.outlet_flow;

    evaluation.residuals = compute_boundary_residuals(
        evaluation.predicted,
        problem.observations);

    // Preserve the observation model's include_initial semantics.
    // At t=0, prescribed boundary values and initial-state
    // variables have a different causal interpretation.
    if (!problem.observations.include_initial) {
        evaluation.residuals.inlet_pressure[0] = 0.0;
        evaluation.residuals.inlet_flow[0] = 0.0;
        evaluation.residuals.outlet_pressure[0] = 0.0;
        evaluation.residuals.outlet_flow[0] = 0.0;
    }

    evaluation.inlet_pressure_cost =
        channel_cost(evaluation.residuals.inlet_pressure);

    evaluation.inlet_flow_cost =
        channel_cost(evaluation.residuals.inlet_flow);

    evaluation.outlet_pressure_cost =
        channel_cost(evaluation.residuals.outlet_pressure);

    evaluation.outlet_flow_cost =
        channel_cost(evaluation.residuals.outlet_flow);

    evaluation.measurement_cost =
        evaluation.inlet_pressure_cost +
        evaluation.inlet_flow_cost +
        evaluation.outlet_pressure_cost +
        evaluation.outlet_flow_cost;

    // No regularizer is enabled in E2.4a.
    evaluation.regularization_cost = 0.0;

    evaluation.total_cost =
        evaluation.measurement_cost +
        evaluation.regularization_cost;

    if (!std::isfinite(evaluation.total_cost))
        throw std::runtime_error("Nonfinite objective value");

    return evaluation;
}

} // namespace tpsf
