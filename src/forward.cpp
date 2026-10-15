#include "tpsf/forward.hpp"
#include <stdexcept>
#include <utility>
namespace tpsf {
Result run_forward(const Configuration& cfg, const InitialState& initial) {
    if (cfg.cells < 2 || cfg.steps < 0)
        throw std::invalid_argument("Invalid cell or step count");
    if (initial.pressure.size() != static_cast<std::size_t>(cfg.cells + 1) ||
        initial.flow.size() != static_cast<std::size_t>(cfg.cells))
        throw std::invalid_argument("Initial-state vector dimension mismatch");
    auto pressure = initial.pressure;
    auto flow = initial.flow;
    auto inlet_flow = initial.inlet_flow;
    Result output;
    const auto count = static_cast<std::size_t>(cfg.steps + 1);
    output.time.resize(count);
    output.inlet_pressure.assign(count, cfg.demand.inlet_pressure);
    output.inlet_flow.resize(count);
    output.outlet_flow.resize(count);
    output.outlet_pressure.resize(count);
    output.linepack.resize(count);
    std::vector<int> iterations(count);
    pipe_sim::IntegrationHistory history{
        output.time, output.inlet_flow, output.outlet_flow,
        output.outlet_pressure, output.linepack, iterations
    };
    pipe_sim::IntegrationWorkspace workspace(cfg.cells);
    const auto result = pipe_sim::integrate_transient(
        cfg.cells, pressure, flow, inlet_flow, cfg.demand, cfg.steps,
        cfg.parameters, cfg.newton, history, workspace);
    output.info = result.info;
    output.completed_steps = result.completed_steps;
    output.max_mass_defect = result.max_mass_defect;
    if (result.info == 0 && result.completed_steps == cfg.steps) {
        output.final_pressure = std::move(pressure);
        output.final_flow = std::move(flow);
        output.final_inlet_flow = inlet_flow;
    }
    return output;
}
} // namespace tpsf
