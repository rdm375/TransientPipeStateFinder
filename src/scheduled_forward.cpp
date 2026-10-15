#include "tpsf/scheduled_forward.hpp"

#include "pipe_sim/transient.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tpsf {

namespace {

std::vector<double> pack_state(
    int n,
    const std::vector<double>& pressure,
    const std::vector<double>& flow,
    double inlet_flow)
{
    std::vector<double> state(2*n + 2);

    state[0] = pressure[0];
    state[1] = inlet_flow;

    for (int i = 0; i < n; ++i) {
        state[2*i + 2] = flow[i];
        state[2*i + 3] = pressure[i + 1];
    }

    return state;
}

Result integrate_scheduled(
    const Configuration& cfg,
    const InitialState& initial,
    const BoundarySchedule& schedule,
    ScoredTrajectory* recorded)
{
    const int n = cfg.cells;
    const int steps = cfg.steps;
    const auto& par = cfg.parameters;

    if (n < 2 || steps < 0)
        throw std::invalid_argument(
            "Invalid integration dimensions");

    const auto count =
        static_cast<std::size_t>(steps) + 1;

    if (initial.pressure.size() !=
            static_cast<std::size_t>(n+1) ||
        initial.flow.size() !=
            static_cast<std::size_t>(n) ||
        schedule.inlet_pressure.size() != count ||
        schedule.outlet_flow.size() != count)
        throw std::invalid_argument(
            "Invalid state or schedule dimensions");

    if (!std::isfinite(par.dt) || par.dt <= 0.0 ||
        !std::isfinite(par.theta) ||
        par.theta < 0.0 || par.theta > 1.0)
        throw std::invalid_argument(
            "Invalid integration parameters");

    for (double value : schedule.inlet_pressure) {
        if (!std::isfinite(value) || value <= 0.0)
            throw std::invalid_argument(
                "Invalid scheduled inlet pressure");
    }

    for (double value : schedule.outlet_flow) {
        if (!std::isfinite(value))
            throw std::invalid_argument(
                "Invalid scheduled outlet flow");
    }

    for (double value : initial.pressure) {
        if (!std::isfinite(value) || value <= 0.0)
            throw std::invalid_argument(
                "Invalid initial pressure");
    }

    for (double value : initial.flow) {
        if (!std::isfinite(value))
            throw std::invalid_argument(
                "Invalid initial flow");
    }

    if (!std::isfinite(initial.inlet_flow))
        throw std::invalid_argument(
            "Invalid initial inlet flow");

    if (initial.pressure[0] !=
        schedule.inlet_pressure[0])
        throw std::invalid_argument(
            "Initial inlet pressure disagrees with boundary schedule");

    auto pressure = initial.pressure;
    auto flow = initial.flow;
    double inlet_flow = initial.inlet_flow;

    pipe_sim::IntegrationWorkspace ws(n);

    Result output;

    output.time.resize(count);
    output.inlet_pressure = schedule.inlet_pressure;
    output.inlet_flow.resize(count);
    output.outlet_flow.resize(count);
    output.outlet_pressure.resize(count);
    output.linepack.resize(count);

    output.time[0] = 0.0;
    output.inlet_flow[0] = inlet_flow;
    output.outlet_flow[0] = schedule.outlet_flow[0];
    output.outlet_pressure[0] = pressure[n];

    output.linepack[0] =
        pipe_sim::linepack_cz(n, pressure, par);

    if (recorded != nullptr) {
        recorded->cells = n;
        recorded->snapshots.clear();
        recorded->snapshots.reserve(count);

        ScoredStateSnapshot snapshot;

        snapshot.step = 0;
        snapshot.time = 0.0;

        snapshot.state =
            pack_state(n, pressure, flow, inlet_flow);

        snapshot.inputs = {
            schedule.inlet_pressure[0],
            schedule.outlet_flow[0]
        };

        snapshot.physical_channels =
            extract_physical_boundary_channels(
                n, snapshot.state, snapshot.inputs);

        snapshot.quality.linepack =
            output.linepack[0];

        snapshot.quality.mass_flow_imbalance =
            inlet_flow - schedule.outlet_flow[0];

        recorded->snapshots.push_back(
            std::move(snapshot));
    }

    double max_mass_defect = 0.0;

    for (int k = 1; k <= steps; ++k) {
        const auto index =
            static_cast<std::size_t>(k);

        const pipe_sim::TransientBoundary bc{
            inlet_flow,
            schedule.outlet_flow[index-1],
            schedule.outlet_flow[index],
            schedule.inlet_pressure[index]
        };

        const auto step = pipe_sim::transient_step(
            n,
            pressure,
            flow,
            bc,
            par,
            cfg.newton,
            ws.state,
            ws.newton);

        if (step.info != 0) {
            output.info = step.info;
            output.completed_steps = k-1;
            output.max_mass_defect = max_mass_defect;
            return output;
        }

        ws.new_pressure[0] = ws.state[0];

        for (int i = 0; i < n; ++i) {
            ws.new_flow[i] = ws.state[2*i+2];
            ws.new_pressure[i+1] = ws.state[2*i+3];
        }

        output.time[index] =
            static_cast<double>(k) * par.dt;

        output.inlet_flow[index] = ws.state[1];
        output.outlet_flow[index] =
            schedule.outlet_flow[index];

        output.outlet_pressure[index] =
            ws.new_pressure[n];

        output.linepack[index] =
            pipe_sim::linepack_cz(
                n, ws.new_pressure, par);

        // Preserve the validated integrator's arithmetic grouping.
        const double balance_error = std::abs(
            (output.linepack[index] -
             output.linepack[index-1]) / par.dt
            - (par.theta *
                (output.inlet_flow[index] -
                 output.outlet_flow[index])
              + (1.0 - par.theta) *
                (output.inlet_flow[index-1] -
                 output.outlet_flow[index-1])));

        max_mass_defect =
            std::max(max_mass_defect, balance_error);

        if (recorded != nullptr) {
            ScoredStateSnapshot snapshot;

            snapshot.step = index;
            snapshot.time = output.time[index];

            snapshot.state = ws.state;

            snapshot.inputs = {
                schedule.inlet_pressure[index],
                schedule.outlet_flow[index]
            };

            snapshot.physical_channels =
                extract_physical_boundary_channels(
                    n, snapshot.state, snapshot.inputs);

            snapshot.quality.linepack =
                output.linepack[index];

            snapshot.quality.mass_flow_imbalance =
                output.inlet_flow[index] -
                output.outlet_flow[index];

            snapshot.quality.mass_conservation_defect =
                balance_error;

            std::vector<double> residual(2*n+2);

            pipe_sim::assemble_residual(
                n,
                snapshot.state,
                pressure,
                flow,
                bc,
                par,
                residual);

            double residual_linf = 0.0;

            for (double value : residual) {
                residual_linf =
                    std::max(residual_linf,
                             std::abs(value));
            }

            snapshot.quality.dynamics_residual_linf =
                residual_linf;

            recorded->snapshots.push_back(
                std::move(snapshot));
        }

        std::copy(
            ws.new_pressure.begin(),
            ws.new_pressure.end(),
            pressure.begin());

        std::copy(
            ws.new_flow.begin(),
            ws.new_flow.end(),
            flow.begin());

        inlet_flow = output.inlet_flow[index];
    }

    output.info = 0;
    output.completed_steps = steps;
    output.max_mass_defect = max_mass_defect;

    output.final_pressure = std::move(pressure);
    output.final_flow = std::move(flow);
    output.final_inlet_flow = inlet_flow;

    return output;
}

} // namespace

Result run_forward_scheduled(
    const Configuration& configuration,
    const InitialState& initial,
    const BoundarySchedule& schedule)
{
    return integrate_scheduled(
        configuration,
        initial,
        schedule,
        nullptr);
}

RecordedForwardResult run_forward_scheduled_recorded(
    const Configuration& configuration,
    const InitialState& initial,
    const BoundarySchedule& schedule)
{
    RecordedForwardResult result;

    result.summary = integrate_scheduled(
        configuration,
        initial,
        schedule,
        &result.trajectory);

    return result;
}

} // namespace tpsf
