#include "tpsf/full_space_constraints.hpp"
#include "tpsf/scheduled_forward.hpp"

#include "pipe_sim/transient_step.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <iostream>
#include <vector>

namespace {

bool close(double a, double b, double tolerance = 1e-10)
{
    return std::abs(a-b) <=
           tolerance*std::max({1.0,std::abs(a),std::abs(b)});
}

void test_observation_independence()
{
    constexpr int n = 2;

    const std::vector<double> state_old{
        5.0e6, 10.0, 10.0, 4.9e6, 10.0, 4.8e6
    };

    const std::vector<double> state_new = state_old;

    const tpsf::PhysicalBoundaryInputs inputs{
        5.0e6, 10.0
    };

    const auto physical =
        tpsf::extract_physical_boundary_channels(
            n, state_new, inputs);

    auto reconstructed = physical;

    const pipe_sim::TransientParameters parameters{
        .dx = 100.0,
        .dt = 1.0,
        .theta = 0.65,
        .diameter = 0.5,
        .area = 0.19634954084936207,
        .temperature = 300.0,
        .z = 0.9,
        .gas_constant = 500.0,
        .viscosity = 1.0e-5,
        .roughness = 1.0e-5
    };

    const auto baseline =
        tpsf::assemble_full_space_step_constraints(
            n,
            state_old,
            state_new,
            inputs,
            inputs,
            reconstructed,
            parameters);

    for (double value : baseline.observation_links)
        assert(close(value, 0.0));

    // Perturb each reconstructed channel independently.
    for (std::size_t channel = 0; channel < 4; ++channel) {
        auto changed = reconstructed;
        changed[channel] += 1.0;

        const auto result =
            tpsf::assemble_full_space_step_constraints(
                n,
                state_old,
                state_new,
                inputs,
                inputs,
                changed,
                parameters);

        for (std::size_t j = 0; j < 4; ++j) {
            assert(close(
                result.observation_links[j],
                j == channel ? 1.0 : 0.0));
        }

        // Observation reconstruction does not modify dynamics.
        assert(result.dynamics == baseline.dynamics);
    }

    assert(baseline.combined().size() ==
           static_cast<std::size_t>(2*n + 6));
}

void test_physical_boundary_linkage()
{
    constexpr int n = 2;

    const std::vector<double> state{
        5.0e6, 10.0, 10.0, 4.9e6, 10.0, 4.8e6
    };

    const tpsf::PhysicalBoundaryInputs inputs{
        5.0e6, 10.0
    };

    const auto reconstructed =
        tpsf::extract_physical_boundary_channels(
            n, state, inputs);

    const pipe_sim::TransientParameters parameters{
        100.0, 1.0, 0.65, 0.5,
        0.19634954084936207, 300.0, 0.9,
        500.0, 1.0e-5, 1.0e-5
    };

    auto changed_inputs = inputs;
    changed_inputs.inlet_pressure += 10.0;
    changed_inputs.outlet_flow += 2.0;

    const auto result =
        tpsf::assemble_full_space_step_constraints(
            n, state, state,
            inputs, changed_inputs,
            reconstructed, parameters);

    assert(close(result.observation_links[0], -10.0));
    assert(close(result.observation_links[3], -2.0));

    // New inlet pressure appears in the solver's boundary equation.
    assert(close(result.dynamics[0], -10.0));

    // New outlet flow appears in the outlet half-cell equation
    // with coefficient theta.
    assert(close(
        result.dynamics.back(),
        0.65*2.0));
}

void test_invalid_dimensions()
{
    constexpr int n = 2;

    const std::vector<double> invalid(3, 0.0);
    const std::vector<double> valid(2*n + 2, 0.0);

    const pipe_sim::TransientParameters parameters{
        100.0, 1.0, 0.65, 0.5,
        0.19634954084936207, 300.0, 0.9,
        500.0, 1.0e-5, 1.0e-5
    };

    bool threw = false;

    try {
        (void)tpsf::assemble_full_space_step_constraints(
            n, invalid, valid,
            {}, {}, {}, parameters);
    } catch (const std::invalid_argument&) {
        threw = true;
    }

    assert(threw);
}


void test_converged_forward_timestep()
{
    constexpr int n = 10;

    const pipe_sim::TransientParameters parameters{
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
    };

    const tpsf::Configuration cfg{
        n,
        1,
        parameters,
        pipe_sim::DemandRamp{
            5.0e6, 10.0, 12.0, 20.0
        },
        pipe_sim::NewtonOptions{
            1.0e-8, 1.0e-10, 30, false
        }
    };

    const tpsf::InitialState initial{
        std::vector<double>(n+1, 5.0e6),
        std::vector<double>(n, 10.0),
        10.0
    };

    tpsf::BoundarySchedule schedule;

    schedule.inlet_pressure = {
        5.0e6,
        5.0e6
    };

    schedule.outlet_flow = {
        10.0,
        10.1
    };

    const auto forward =
        tpsf::run_forward_scheduled(
            cfg, initial, schedule);

    if (forward.info != 0 ||
        forward.completed_steps != 1) {
        throw std::runtime_error(
            "Converged-timestep fixture failed");
    }

    // Construct the old state using the solver's
    // interleaved physical-state ordering.
    std::vector<double> state_old(2*n+2);

    state_old[0] = initial.pressure[0];
    state_old[1] = initial.inlet_flow;

    for (int i = 0; i < n; ++i) {
        state_old[2*i+2] = initial.flow[i];
        state_old[2*i+3] = initial.pressure[i+1];
    }

    // The one-step integration result is x_1.
    std::vector<double> state_new(2*n+2);

    state_new[0] = forward.final_pressure[0];
    state_new[1] = forward.final_inlet_flow;

    for (int i = 0; i < n; ++i) {
        state_new[2*i+2] = forward.final_flow[i];
        state_new[2*i+3] =
            forward.final_pressure[i+1];
    }

    const tpsf::PhysicalBoundaryInputs inputs_old{
        schedule.inlet_pressure[0],
        schedule.outlet_flow[0]
    };

    const tpsf::PhysicalBoundaryInputs inputs_new{
        schedule.inlet_pressure[1],
        schedule.outlet_flow[1]
    };

    const auto reconstructed =
        tpsf::extract_physical_boundary_channels(
            n, state_new, inputs_new);

    const auto constraints =
        tpsf::assemble_full_space_step_constraints(
            n,
            state_old,
            state_new,
            inputs_old,
            inputs_new,
            reconstructed,
            parameters);

    // Independent direct assembly using the solver API.
    const pipe_sim::TransientBoundary bc{
        initial.inlet_flow,
        schedule.outlet_flow[0],
        schedule.outlet_flow[1],
        schedule.inlet_pressure[1]
    };

    std::vector<double> direct_residual(2*n+2);

    pipe_sim::assemble_residual(
        n,
        state_new,
        initial.pressure,
        initial.flow,
        bc,
        parameters,
        direct_residual);

    double maximum_residual = 0.0;
    double maximum_difference = 0.0;
    double maximum_link = 0.0;

    for (std::size_t i = 0;
         i < constraints.dynamics.size(); ++i) {

        maximum_residual = std::max(
            maximum_residual,
            std::abs(constraints.dynamics[i]));

        maximum_difference = std::max(
            maximum_difference,
            std::abs(
                constraints.dynamics[i] -
                direct_residual[i]));
    }

    for (double value : constraints.observation_links) {
        maximum_link = std::max(
            maximum_link,
            std::abs(value));
    }

    std::cout
        << "Converged timestep:\n"
        << "  max dynamics residual = "
        << maximum_residual << '\n'
        << "  max assembly difference = "
        << maximum_difference << '\n'
        << "  max observation link = "
        << maximum_link << '\n';

    if (!std::isfinite(maximum_residual) ||
        !std::isfinite(maximum_difference) ||
        !std::isfinite(maximum_link)) {
        throw std::runtime_error(
            "Nonfinite converged-timestep residual");
    }

    // Fixture-specific regression threshold.
    // Measured reference residual: approximately 4.16e-12.
    if (maximum_residual >= 1.0e-9) {
        throw std::runtime_error(
            "Converged dynamics residual exceeds tolerance");
    }

    if (maximum_difference > 1.0e-12) {
        throw std::runtime_error(
            "Full-space dynamics differ from solver residual");
    }

    if (maximum_link != 0.0) {
        throw std::runtime_error(
            "Physical observation links are inconsistent");
    }

    if (constraints.combined().size() !=
        static_cast<std::size_t>(2*n+6)) {
        throw std::runtime_error(
            "Incorrect combined constraint dimensions");
    }
}

} // namespace

int main()
{
    test_observation_independence();
    test_physical_boundary_linkage();
    test_invalid_dimensions();
    test_converged_forward_timestep();
}
