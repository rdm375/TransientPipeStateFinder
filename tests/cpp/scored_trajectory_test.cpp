#include "tpsf/scheduled_forward.hpp"
#include "tpsf/full_space_trajectory.hpp"

#include "pipe_sim/transient.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(
    double actual,
    double expected,
    double atol,
    double rtol,
    const std::string& message)
{
    require(
        std::isfinite(actual) &&
        std::isfinite(expected),
        message + ": nonfinite value");

    require(
        std::abs(actual - expected) <=
            atol + rtol *
                std::max(std::abs(actual),
                         std::abs(expected)),
        message + ": numerical disagreement");
}

tpsf::Configuration configuration(int steps = 12)
{
    constexpr int n = 10;

    return {
        n,
        steps,
        pipe_sim::TransientParameters{
            100.0,
            1.0,
            0.65,
            0.5,
            0.19634954084936207,
            288.15,
            0.9,
            518.3,
            1.1e-5,
            1.0e-5
        },
        pipe_sim::DemandRamp{
            5.0e6,
            10.0,
            12.0,
            20.0
        },
        pipe_sim::NewtonOptions{
            1.0e-8,
            1.0e-10,
            30,
            false
        }
    };
}

tpsf::InitialState initial_state(int n)
{
    return {
        std::vector<double>(n+1, 5.0e6),
        std::vector<double>(n, 10.0),
        10.0
    };
}

tpsf::BoundarySchedule schedule_for(
    const tpsf::Configuration& cfg)
{
    tpsf::BoundarySchedule schedule;

    for (int k = 0; k <= cfg.steps; ++k) {
        const double time =
            k * cfg.parameters.dt;

        schedule.inlet_pressure.push_back(
            cfg.demand.inlet_pressure);

        schedule.outlet_flow.push_back(
            pipe_sim::demand_ramp(
                time,
                cfg.demand.initial_outlet_flow,
                cfg.demand.final_outlet_flow,
                cfg.demand.ramp_duration));
    }

    return schedule;
}

void test_recorded_forward_parity()
{
    const auto cfg = configuration();
    const auto initial = initial_state(cfg.cells);
    const auto schedule = schedule_for(cfg);

    const auto reference =
        tpsf::run_forward_scheduled(
            cfg, initial, schedule);

    const auto recorded =
        tpsf::run_forward_scheduled_recorded(
            cfg, initial, schedule);

    require(reference.info == 0,
            "Reference integration failed");

    require(recorded.summary.info == 0,
            "Recorded integration failed");

    require(
        recorded.summary.completed_steps == cfg.steps,
        "Recorded integration incomplete");

    require(
        recorded.trajectory.cells == cfg.cells,
        "Recorded cell count incorrect");

    require(
        recorded.trajectory.snapshots.size() ==
            static_cast<std::size_t>(cfg.steps+1),
        "Incorrect snapshot count");

    const auto& summary = recorded.summary;

    require(
        summary.time == reference.time,
        "Time history changed");

    require(
        summary.inlet_pressure ==
            reference.inlet_pressure,
        "Inlet-pressure history changed");

    require(
        summary.inlet_flow ==
            reference.inlet_flow,
        "Inlet-flow history changed");

    require(
        summary.outlet_flow ==
            reference.outlet_flow,
        "Outlet-flow history changed");

    require(
        summary.outlet_pressure ==
            reference.outlet_pressure,
        "Outlet-pressure history changed");

    require(
        summary.linepack ==
            reference.linepack,
        "Linepack history changed");

    require(
        summary.final_pressure ==
            reference.final_pressure,
        "Final pressure changed");

    require(
        summary.final_flow ==
            reference.final_flow,
        "Final flow changed");

    near(
        summary.final_inlet_flow,
        reference.final_inlet_flow,
        0.0,
        0.0,
        "Final inlet flow changed");

    near(
        summary.max_mass_defect,
        reference.max_mass_defect,
        0.0,
        0.0,
        "Maximum conservation defect changed");

    std::cout
        << "PASS recorded forward parity\n";
}

void test_snapshot_contents()
{
    const auto cfg = configuration();
    const auto initial = initial_state(cfg.cells);
    const auto schedule = schedule_for(cfg);

    const auto result =
        tpsf::run_forward_scheduled_recorded(
            cfg, initial, schedule);

    require(result.summary.info == 0,
            "Integration failed");

    const auto& snapshots =
        result.trajectory.snapshots;

    const int n = cfg.cells;

    for (int k = 0; k <= cfg.steps; ++k) {
        const auto& snapshot = snapshots[k];

        require(
            snapshot.step ==
                static_cast<std::size_t>(k),
            "Snapshot index incorrect");

        near(
            snapshot.time,
            k * cfg.parameters.dt,
            0.0,
            0.0,
            "Snapshot time incorrect");

        require(
            snapshot.state.size() ==
                static_cast<std::size_t>(2*n+2),
            "Snapshot state dimension incorrect");

        near(
            snapshot.inputs.inlet_pressure,
            schedule.inlet_pressure[k],
            0.0,
            0.0,
            "Snapshot inlet pressure incorrect");

        near(
            snapshot.inputs.outlet_flow,
            schedule.outlet_flow[k],
            0.0,
            0.0,
            "Snapshot outlet flow incorrect");

        const auto physical =
            tpsf::extract_physical_boundary_channels(
                n,
                snapshot.state,
                snapshot.inputs);

        for (std::size_t c = 0; c < 4; ++c) {
            near(
                snapshot.physical_channels[c],
                physical[c],
                0.0,
                0.0,
                "Physical boundary channel incorrect");
        }

        near(
            snapshot.physical_channels[0],
            result.summary.inlet_pressure[k],
            0.0,
            0.0,
            "Recorded inlet-pressure history mismatch");

        near(
            snapshot.physical_channels[1],
            result.summary.inlet_flow[k],
            0.0,
            0.0,
            "Recorded inlet-flow history mismatch");

        near(
            snapshot.physical_channels[2],
            result.summary.outlet_pressure[k],
            0.0,
            0.0,
            "Recorded outlet-pressure history mismatch");

        near(
            snapshot.physical_channels[3],
            result.summary.outlet_flow[k],
            0.0,
            0.0,
            "Recorded outlet-flow history mismatch");

        require(
            snapshot.quality.linepack.has_value(),
            "Linepack score missing");

        near(
            *snapshot.quality.linepack,
            result.summary.linepack[k],
            0.0,
            0.0,
            "Recorded linepack mismatch");

        require(
            snapshot.quality.mass_flow_imbalance.has_value(),
            "Physical imbalance score missing");

        near(
            *snapshot.quality.mass_flow_imbalance,
            result.summary.inlet_flow[k] -
                result.summary.outlet_flow[k],
            0.0,
            0.0,
            "Physical imbalance incorrect");

        require(
            !snapshot.quality.measurement_weighted_rms,
            "Measurement score unexpectedly available");

        require(
            !snapshot.quality.measurement_max_normalized_error,
            "Maximum measurement error unexpectedly available");

        require(
            snapshot.quality.measurement_count == 0,
            "Unexpected measurement count");

        require(
            !snapshot.quality.state_uncertainty_rms,
            "State uncertainty unexpectedly available");

        require(
            !snapshot.quality.state_uncertainty_max,
            "Maximum uncertainty unexpectedly available");
    }

    const auto& first = snapshots.front();

    near(
        first.state[0],
        initial.pressure[0],
        0.0,
        0.0,
        "Initial pressure incorrect");

    near(
        first.state[1],
        initial.inlet_flow,
        0.0,
        0.0,
        "Initial inlet flow incorrect");

    for (int i = 0; i < n; ++i) {
        near(
            first.state[2*i+2],
            initial.flow[i],
            0.0,
            0.0,
            "Initial interior flow incorrect");

        near(
            first.state[2*i+3],
            initial.pressure[i+1],
            0.0,
            0.0,
            "Initial interior pressure incorrect");
    }

    require(
        !first.quality.dynamics_residual_linf,
        "Initial dynamics residual should be unavailable");

    require(
        !first.quality.mass_conservation_defect,
        "Initial conservation defect should be unavailable");

    const auto& last = snapshots.back();

    near(
        last.state[1],
        result.summary.final_inlet_flow,
        0.0,
        0.0,
        "Final recorded inlet flow incorrect");

    for (int i = 0; i < n; ++i) {
        near(
            last.state[2*i+2],
            result.summary.final_flow[i],
            0.0,
            0.0,
            "Final recorded interior flow incorrect");

        near(
            last.state[2*i+3],
            result.summary.final_pressure[i+1],
            0.0,
            0.0,
            "Final recorded pressure incorrect");
    }

    std::cout
        << "PASS snapshot contents\n";
}

void test_numerical_scores()
{
    const auto cfg = configuration();
    const auto initial = initial_state(cfg.cells);
    const auto schedule = schedule_for(cfg);

    const auto result =
        tpsf::run_forward_scheduled_recorded(
            cfg, initial, schedule);

    require(result.summary.info == 0,
            "Integration failed");

    const auto& snapshots =
        result.trajectory.snapshots;

    double maximum_defect = 0.0;

    for (int k = 1; k <= cfg.steps; ++k) {
        const auto& previous = snapshots[k-1];
        const auto& current = snapshots[k];

        require(
            current.quality.mass_conservation_defect.has_value(),
            "Missing conservation defect");

        require(
            current.quality.dynamics_residual_linf.has_value(),
            "Missing dynamics residual");

        const double imbalance_new =
            result.summary.inlet_flow[k] -
            result.summary.outlet_flow[k];

        const double imbalance_old =
            result.summary.inlet_flow[k-1] -
            result.summary.outlet_flow[k-1];

        const double defect = std::abs(
            (result.summary.linepack[k] -
             result.summary.linepack[k-1]) /
                cfg.parameters.dt
            - (cfg.parameters.theta *
                imbalance_new
              + (1.0 - cfg.parameters.theta) *
                imbalance_old));

        near(
            *current.quality.mass_conservation_defect,
            defect,
            0.0,
            0.0,
            "Conservation defect mismatch");

        maximum_defect =
            std::max(maximum_defect, defect);

        std::vector<double> old_pressure(
            cfg.cells+1);

        std::vector<double> old_flow(
            cfg.cells);

        old_pressure[0] = previous.state[0];

        for (int i = 0; i < cfg.cells; ++i) {
            old_flow[i] =
                previous.state[2*i+2];

            old_pressure[i+1] =
                previous.state[2*i+3];
        }

        const pipe_sim::TransientBoundary boundary{
            previous.state[1],
            previous.inputs.outlet_flow,
            current.inputs.outlet_flow,
            current.inputs.inlet_pressure
        };

        std::vector<double> residual(
            2*cfg.cells+2);

        pipe_sim::assemble_residual(
            cfg.cells,
            current.state,
            old_pressure,
            old_flow,
            boundary,
            cfg.parameters,
            residual);

        double residual_linf = 0.0;

        for (double value : residual) {
            residual_linf =
                std::max(
                    residual_linf,
                    std::abs(value));
        }

        near(
            *current.quality.dynamics_residual_linf,
            residual_linf,
            0.0,
            0.0,
            "Dynamics residual score mismatch");

        require(
            std::isfinite(residual_linf),
            "Nonfinite dynamics residual");

        require(
            residual_linf < 1.0e-5,
            "Converged dynamics residual unexpectedly large");
    }

    near(
        maximum_defect,
        result.summary.max_mass_defect,
        0.0,
        0.0,
        "Maximum mass defect mismatch");

    std::cout
        << "PASS numerical quality scores\n";
}

void test_full_space_consistency()
{
    const auto cfg = configuration();
    const auto initial = initial_state(cfg.cells);
    const auto schedule = schedule_for(cfg);

    const auto recorded =
        tpsf::run_forward_scheduled_recorded(
            cfg, initial, schedule);

    require(recorded.summary.info == 0,
            "Integration failed");

    tpsf::FullSpaceTrajectory trajectory;
    trajectory.n = cfg.cells;

    for (const auto& snapshot :
         recorded.trajectory.snapshots) {

        trajectory.states.push_back(
            snapshot.state);

        trajectory.inputs.push_back(
            snapshot.inputs);

        trajectory.reconstructed.push_back(
            snapshot.physical_channels);
    }

    const auto constraints =
        tpsf::assemble_full_space_trajectory_constraints(
            trajectory,
            cfg.parameters);

    require(
        constraints.dynamics.size() ==
            static_cast<std::size_t>(cfg.steps),
        "Full-space dynamics count incorrect");

    require(
        constraints.observation_links.size() ==
            static_cast<std::size_t>(cfg.steps+1),
        "Full-space link count incorrect");

    for (int k = 0; k < cfg.steps; ++k) {
        double residual_linf = 0.0;

        for (double value : constraints.dynamics[k]) {
            residual_linf =
                std::max(
                    residual_linf,
                    std::abs(value));
        }

        near(
            residual_linf,
            *recorded.trajectory.snapshots[k+1]
                 .quality.dynamics_residual_linf,
            0.0,
            0.0,
            "Full-space dynamics score mismatch");
    }

    for (const auto& link :
         constraints.observation_links) {

        for (double value : link) {
            near(
                value,
                0.0,
                0.0,
                0.0,
                "Physical observation link nonzero");
        }
    }

    std::cout
        << "PASS full-space trajectory consistency\n";
}

void test_zero_steps()
{
    const auto cfg = configuration(0);
    const auto initial = initial_state(cfg.cells);
    const auto schedule = schedule_for(cfg);

    const auto result =
        tpsf::run_forward_scheduled_recorded(
            cfg, initial, schedule);

    require(
        result.summary.info == 0,
        "Zero-step integration failed");

    require(
        result.summary.completed_steps == 0,
        "Zero-step completed count incorrect");

    require(
        result.trajectory.snapshots.size() == 1,
        "Zero-step trajectory must contain initial state");

    require(
        !result.trajectory.snapshots[0]
             .quality.dynamics_residual_linf,
        "Zero-step dynamics score should be unavailable");

    require(
        !result.trajectory.snapshots[0]
             .quality.mass_conservation_defect,
        "Zero-step conservation score should be unavailable");

    std::cout
        << "PASS zero-step trajectory\n";
}

} // namespace

int main()
{
    try {
        test_recorded_forward_parity();
        test_snapshot_contents();
        test_numerical_scores();
        test_full_space_consistency();
        test_zero_steps();

        std::cout
            << "PASS E2.4b.7 scored trajectory recording\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E2.4b.7: "
            << error.what()
            << '\n';

        return 1;
    }
}
