#pragma once

#include "tpsf/forward.hpp"

#include <vector>

namespace tpsf {

struct BoundarySchedule {
    // Samples at t_k = k * parameters.dt, k = 0,...,steps.
    std::vector<double> inlet_pressure;
    std::vector<double> outlet_flow;
};

Result run_forward_scheduled(
    const Configuration& configuration,
    const InitialState& initial,
    const BoundarySchedule& schedule);

} // namespace tpsf
