#pragma once

#include "tpsf/forward.hpp"
#include "tpsf/scored_trajectory.hpp"

#include <vector>

namespace tpsf {

struct BoundarySchedule {
    std::vector<double> inlet_pressure;
    std::vector<double> outlet_flow;
};

struct RecordedForwardResult {
    Result summary;
    ScoredTrajectory trajectory;
};

Result run_forward_scheduled(
    const Configuration& configuration,
    const InitialState& initial,
    const BoundarySchedule& schedule);

RecordedForwardResult run_forward_scheduled_recorded(
    const Configuration& configuration,
    const InitialState& initial,
    const BoundarySchedule& schedule);

} // namespace tpsf
