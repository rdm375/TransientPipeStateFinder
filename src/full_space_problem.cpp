#include "tpsf/full_space_problem.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace tpsf {
namespace {

std::size_t state_width(int n)
{
    if (n < 2)
        throw std::invalid_argument("At least two cells required");

    return static_cast<std::size_t>(2*n+2);
}

std::size_t block_width(int n)
{
    return state_width(n)+6;
}

void validate_trajectory(const FullSpaceTrajectory& t)
{
    const auto width = state_width(t.n);
    const auto count = t.states.size();

    if (count < 2 ||
        t.inputs.size() != count ||
        t.reconstructed.size() != count)
        throw std::invalid_argument(
            "Invalid trajectory time dimensions");

    for (const auto& state : t.states)
        if (state.size() != width)
            throw std::invalid_argument(
                "Invalid trajectory state width");
}

std::vector<double> flatten(const FullSpaceTrajectory& t)
{
    validate_trajectory(t);

    std::vector<double> flat;
    flat.reserve(t.states.size()*block_width(t.n));

    for (std::size_t k = 0; k < t.states.size(); ++k) {
        flat.insert(flat.end(),
                    t.states[k].begin(),
                    t.states[k].end());

        flat.push_back(t.inputs[k].inlet_pressure);
        flat.push_back(t.inputs[k].outlet_flow);

        flat.insert(flat.end(),
                    t.reconstructed[k].begin(),
                    t.reconstructed[k].end());
    }

    return flat;
}

FullSpaceTrajectory unflatten(
    const std::vector<double>& flat,
    int n,
    std::size_t count)
{
    const auto width = state_width(n);

    if (flat.size() != count*(width+6))
        throw std::invalid_argument(
            "Incorrect flattened trajectory size");

    FullSpaceTrajectory result;
    result.n = n;

    std::size_t offset = 0;

    for (std::size_t k = 0; k < count; ++k) {
        result.states.emplace_back(
            flat.begin()+offset,
            flat.begin()+offset+width);

        offset += width;

        PhysicalBoundaryInputs inputs;
        inputs.inlet_pressure = flat[offset++];
        inputs.outlet_flow = flat[offset++];
        result.inputs.push_back(inputs);

        BoundaryVector reconstructed{};

        for (std::size_t c = 0; c < 4; ++c)
            reconstructed[c] = flat[offset++];

        result.reconstructed.push_back(reconstructed);
    }

    return result;
}

std::vector<double> flatten(
    const FullSpaceTrajectoryAdjoint& a,
    int n)
{
    FullSpaceTrajectory temporary;
    temporary.n = n;
    temporary.states = a.states;
    temporary.inputs = a.inputs;
    temporary.reconstructed = a.reconstructed;

    return flatten(temporary);
}

FullSpaceTrajectoryDirection direction_from_flat(
    const std::vector<double>& flat,
    int n,
    std::size_t count)
{
    const auto t = unflatten(flat, n, count);

    FullSpaceTrajectoryDirection result;
    result.states = t.states;
    result.inputs = t.inputs;
    result.reconstructed = t.reconstructed;

    return result;
}

FullSpaceTrajectoryConstraints dual_from_flat(
    const std::vector<double>& flat,
    int n,
    std::size_t count)
{
    const auto width = state_width(n);
    const auto steps = count-1;

    if (flat.size() != steps*width+4*count)
        throw std::invalid_argument(
            "Incorrect constraint dual dimension");

    FullSpaceTrajectoryConstraints result;

    std::size_t offset = 0;

    for (std::size_t k = 0; k < steps; ++k) {
        result.dynamics.emplace_back(
            flat.begin()+offset,
            flat.begin()+offset+width);

        offset += width;
    }

    for (std::size_t k = 0; k < count; ++k) {
        BoundaryVector block{};

        for (std::size_t c = 0; c < 4; ++c)
            block[c] = flat[offset++];

        result.observation_links.push_back(block);
    }

    return result;
}

} // namespace

FullSpaceProblem::FullSpaceProblem(
    FullSpaceProblemConfiguration configuration)
    : configuration_(std::move(configuration))
{
    validate_trajectory(configuration_.reference);

    reference_flat_ = flatten(configuration_.reference);

    if (configuration_.variables.size() !=
        reference_flat_.size())
        throw std::invalid_argument(
            "Variable configuration size mismatch");

    if (configuration_.measurements.size() !=
        configuration_.reference.states.size())
        throw std::invalid_argument(
            "Measurement time count mismatch");

    const auto count = configuration_.reference.states.size();
    const auto width = state_width(configuration_.reference.n);

    const auto expected_constraints =
        (count-1)*width+4*count;

    if (configuration_.constraint_scales.size() !=
        expected_constraints)
        throw std::invalid_argument(
            "Constraint scale count mismatch");

    for (double value : reference_flat_)
        if (!std::isfinite(value))
            throw std::invalid_argument(
                "Nonfinite reference trajectory");

    for (std::size_t i = 0;
         i < configuration_.variables.size(); ++i) {
        const auto& variable = configuration_.variables[i];

        if (!std::isfinite(variable.scale) ||
            !(variable.scale > 0.0))
            throw std::invalid_argument(
                "Invalid variable scale");

        if (variable.free)
            free_indices_.push_back(i);
    }

    for (double scale : configuration_.constraint_scales)
        if (!std::isfinite(scale) || !(scale > 0.0))
            throw std::invalid_argument(
                "Invalid constraint scale");

    // Validate the measurement dataset, including the
    // possibility that no channels are available.
    evaluate_full_space_measurement_objective(
        configuration_.reference.reconstructed,
        configuration_.measurements);
}

std::size_t FullSpaceProblem::physical_variable_count() const noexcept
{
    return reference_flat_.size();
}

std::size_t FullSpaceProblem::free_variable_count() const noexcept
{
    return free_indices_.size();
}

std::size_t FullSpaceProblem::constraint_count() const noexcept
{
    return configuration_.constraint_scales.size();
}

FullSpaceTrajectory FullSpaceProblem::trajectory(
    const std::vector<double>& u) const
{
    if (u.size() != free_indices_.size())
        throw std::invalid_argument(
            "Incorrect optimization coordinate count");

    auto flat = reference_flat_;

    for (std::size_t j = 0; j < u.size(); ++j) {
        if (!std::isfinite(u[j]))
            throw std::invalid_argument(
                "Nonfinite optimization coordinate");

        const auto i = free_indices_[j];

        flat[i] += configuration_.variables[i].scale*u[j];

        if (!std::isfinite(flat[i]))
            throw std::invalid_argument(
                "Nonfinite transformed physical variable");
    }

    return unflatten(
        flat,
        configuration_.reference.n,
        configuration_.reference.states.size());
}

std::vector<double> FullSpaceProblem::coordinates(
    const FullSpaceTrajectory& t) const
{
    if (t.n != configuration_.reference.n ||
        t.states.size() != configuration_.reference.states.size())
        throw std::invalid_argument(
            "Incompatible trajectory geometry");

    const auto flat = flatten(t);

    std::vector<double> u;
    u.reserve(free_indices_.size());

    std::size_t free_position = 0;

    for (std::size_t i = 0; i < flat.size(); ++i) {
        if (!std::isfinite(flat[i]))
            throw std::invalid_argument(
                "Nonfinite physical trajectory");

        if (free_position < free_indices_.size() &&
            free_indices_[free_position] == i) {
            u.push_back(
                (flat[i]-reference_flat_[i]) /
                configuration_.variables[i].scale);

            if (!std::isfinite(u.back()))
                throw std::invalid_argument(
                    "Nonfinite recovered coordinate");

            ++free_position;
        } else if (flat[i] != reference_flat_[i]) {
            throw std::invalid_argument(
                "Fixed variable differs from reference");
        }
    }

    return u;
}

std::vector<double> FullSpaceProblem::expand_direction(
    const std::vector<double>& direction) const
{
    if (direction.size() != free_indices_.size())
        throw std::invalid_argument(
            "Incorrect direction dimension");

    std::vector<double> flat(reference_flat_.size(), 0.0);

    for (std::size_t j = 0; j < direction.size(); ++j) {
        if (!std::isfinite(direction[j]))
            throw std::invalid_argument(
                "Nonfinite optimization direction");

        const auto i = free_indices_[j];

        flat[i] =
            configuration_.variables[i].scale*direction[j];

        if (!std::isfinite(flat[i]))
            throw std::invalid_argument(
                "Nonfinite physical direction");
    }

    return flat;
}

std::vector<double> FullSpaceProblem::project_gradient(
    const std::vector<double>& physical_gradient) const
{
    if (physical_gradient.size() != reference_flat_.size())
        throw std::invalid_argument(
            "Incorrect physical gradient dimension");

    std::vector<double> result;
    result.reserve(free_indices_.size());

    for (const auto i : free_indices_) {
        const double value =
            configuration_.variables[i].scale *
            physical_gradient[i];

        if (!std::isfinite(value))
            throw std::invalid_argument(
                "Nonfinite scaled gradient");

        result.push_back(value);
    }

    return result;
}

FullSpaceProblemEvaluation FullSpaceProblem::evaluate(
    const std::vector<double>& u) const
{
    const auto t = trajectory(u);

    FullSpaceProblemEvaluation result;

    result.measurement =
        evaluate_full_space_measurement_objective(
            t.reconstructed,
            configuration_.measurements);

    FullSpaceTrajectoryAdjoint gradient;

    const auto count = t.states.size();
    const auto width = state_width(t.n);

    gradient.states.assign(
        count, std::vector<double>(width, 0.0));

    gradient.inputs.resize(count);
    gradient.reconstructed = result.measurement.gradient;

    result.objective_gradient =
        project_gradient(flatten(gradient, t.n));

    result.constraints =
        assemble_full_space_trajectory_constraints(
            t, configuration_.parameters).combined();

    for (std::size_t i = 0;
         i < result.constraints.size(); ++i) {
        result.constraints[i] *=
            configuration_.constraint_scales[i];
    }

    return result;
}

std::vector<double> FullSpaceProblem::constraint_jvp(
    const std::vector<double>& u,
    const std::vector<double>& direction) const
{
    const auto t = trajectory(u);

    const auto physical_direction =
        expand_direction(direction);

    const auto structured_direction =
        direction_from_flat(
            physical_direction,
            t.n,
            t.states.size());

    auto result =
        apply_full_space_trajectory_jvp(
            t,
            structured_direction,
            configuration_.parameters).combined();

    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] *= configuration_.constraint_scales[i];

    return result;
}

std::vector<double> FullSpaceProblem::constraint_vjp(
    const std::vector<double>& u,
    const std::vector<double>& dual) const
{
    const auto t = trajectory(u);

    if (dual.size() != constraint_count())
        throw std::invalid_argument(
            "Incorrect constraint dual dimension");

    std::vector<double> physical_dual(dual.size());

    for (std::size_t i = 0; i < dual.size(); ++i) {
        if (!std::isfinite(dual[i]))
            throw std::invalid_argument(
                "Nonfinite constraint dual");

        physical_dual[i] =
            configuration_.constraint_scales[i]*dual[i];

        if (!std::isfinite(physical_dual[i]))
            throw std::invalid_argument(
                "Nonfinite scaled constraint dual");
    }

    const auto structured_dual =
        dual_from_flat(
            physical_dual,
            t.n,
            t.states.size());

    const auto physical_gradient =
        apply_full_space_trajectory_vjp(
            t,
            structured_dual,
            configuration_.parameters);

    return project_gradient(
        flatten(physical_gradient, t.n));
}

} // namespace tpsf
