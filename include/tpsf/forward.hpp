#pragma once
#include "pipe_sim/integrate.hpp"
#include <vector>
namespace tpsf {
struct InitialState {
    std::vector<double> pressure;
    std::vector<double> flow;
    double inlet_flow{};
};
struct Configuration {
    int cells{};
    int steps{};
    pipe_sim::TransientParameters parameters;
    pipe_sim::DemandRamp demand;
    pipe_sim::NewtonOptions newton;
};
struct Result {
    int info{};
    int completed_steps{};
    double max_mass_defect{};
    std::vector<double> time, inlet_pressure, inlet_flow;
    std::vector<double> outlet_flow, outlet_pressure, linepack;
    std::vector<double> final_pressure, final_flow;
    double final_inlet_flow{};
};
Result run_forward(const Configuration&, const InitialState&);
} // namespace tpsf
