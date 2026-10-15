#pragma once

#include "tpsf/boundary_residuals.hpp"
#include "tpsf/scaled_initial_state.hpp"
#include "tpsf/scheduled_forward.hpp"

#include <cstddef>
#include <vector>

namespace tpsf {

// Layout:
// [scaled initial state: 2*cells+1]
// [scaled inlet-pressure corrections: steps+1]
// [scaled outlet-flow corrections: steps+1]
struct ObjectiveLayout {
    std::size_t initial_size{};
    std::size_t boundary_size{};
    std::size_t inlet_pressure_offset{};
    std::size_t outlet_flow_offset{};
    std::size_t total_size{};
};

ObjectiveLayout objective_layout(
    const Configuration& configuration);

struct ObjectiveProblem {
    Configuration configuration;
    InitialState baseline_initial;
    BoundaryObservationRecord observations;
    InitialStateScales state_scales{};
};

struct ObjectiveEvaluation {
    double measurement_cost{};
    double regularization_cost{};
    double total_cost{};

    double inlet_pressure_cost{};
    double inlet_flow_cost{};
    double outlet_pressure_cost{};
    double outlet_flow_cost{};

    BoundaryResiduals residuals;
    BoundaryMeasurements predicted;
    Result forward;
};

std::vector<double> zero_objective_parameters(
    const ObjectiveProblem& problem);

ObjectiveEvaluation evaluate_objective(
    const ObjectiveProblem& problem,
    const std::vector<double>& parameters);

} // namespace tpsf
