#pragma once

#include "tpsf/full_space_trajectory_derivatives.hpp"
#include "tpsf/full_space_measurement_objective.hpp"

#include <cstddef>
#include <vector>

namespace tpsf {

// Configuration of one physical optimization variable.
//
// If free == false, the variable is fixed at its reference
// value and does not receive an optimization coordinate.
struct ProblemVariable {
    bool free = true;
    double scale = 1.0;
};

struct FullSpaceProblemConfiguration {
    FullSpaceTrajectory reference;

    // One entry per physical variable in the ordering:
    // [state, inlet pressure, outlet flow, reconstructed]
    // at each time level.
    std::vector<ProblemVariable> variables;

    // One positive multiplier per equality constraint.
    // Ordering matches FullSpaceTrajectoryConstraints::combined().
    std::vector<double> constraint_scales;

    std::vector<BoundaryMeasurement> measurements;

    pipe_sim::TransientParameters parameters;
};

struct FullSpaceProblemEvaluation {
    MeasurementObjectiveResult measurement;
    std::vector<double> objective_gradient;
    std::vector<double> constraints;
};

class FullSpaceProblem {
public:
    explicit FullSpaceProblem(
        FullSpaceProblemConfiguration configuration);

    std::size_t physical_variable_count() const noexcept;
    std::size_t free_variable_count() const noexcept;
    std::size_t constraint_count() const noexcept;

    // Reference trajectory corresponds to u = 0.
    FullSpaceTrajectory trajectory(
        const std::vector<double>& u) const;

    // Inverse transformation, for free coordinates.
    // Fixed coordinates must match the reference exactly.
    std::vector<double> coordinates(
        const FullSpaceTrajectory& trajectory) const;

    FullSpaceProblemEvaluation evaluate(
        const std::vector<double>& u) const;

    // Normalized residuals for available measurements only.
    // Ordering: time level, then boundary channel.
    std::vector<double> measurement_residual(
        const std::vector<double>& u) const;

    // Exact normalized measurement-residual JVP.
    std::vector<double> measurement_jvp(
        const std::vector<double>& u,
        const std::vector<double>& direction) const;

    std::vector<double> constraint_jvp(
        const std::vector<double>& u,
        const std::vector<double>& direction) const;

    std::vector<double> constraint_vjp(
        const std::vector<double>& u,
        const std::vector<double>& dual) const;

private:
    FullSpaceProblemConfiguration configuration_;
    std::vector<double> reference_flat_;
    std::vector<std::size_t> free_indices_;

    std::vector<double> expand_direction(
        const std::vector<double>& direction) const;

    std::vector<double> project_gradient(
        const std::vector<double>& physical_gradient) const;
};

} // namespace tpsf
