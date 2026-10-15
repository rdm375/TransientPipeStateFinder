#include "tpsf/forecast_attenuation.hpp"
#include "tpsf/scheduled_forward.hpp"
#include "pipe_sim/physics.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace tpsf;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

std::vector<double> as_std(const Eigen::VectorXd& x)
{
    return std::vector<double>(x.data(), x.data()+x.size());
}

void near(
    double actual,
    double expected,
    double atol,
    double rtol,
    const std::string& message)
{
    require(
        std::isfinite(actual) && std::isfinite(expected),
        message + ": nonfinite");

    require(
        std::abs(actual-expected) <=
            atol + rtol*std::max(
                std::abs(actual), std::abs(expected)),
        message + ": mismatch");
}

FullSpaceTrajectory recorded_trajectory(
    const RecordedForwardResult& recorded)
{
    FullSpaceTrajectory result;
    result.n = recorded.trajectory.cells;

    for (const auto& snapshot :
         recorded.trajectory.snapshots) {
        result.states.push_back(snapshot.state);
        result.inputs.push_back(snapshot.inputs);
        result.reconstructed.push_back(
            snapshot.physical_channels);
    }

    return result;
}

void test_direction_mapping()
{
    constexpr int n = 2;

    Configuration configuration{
        n,
        1,
        pipe_sim::TransientParameters{
            100.0, 1.0, 0.65,
            0.5, 0.19634954084936207,
            288.15, 0.9, 518.3,
            1.1e-5, 1.0e-5
        },
        pipe_sim::DemandRamp{
            5.0e6, 10.0, 10.1, 1.0
        },
        pipe_sim::NewtonOptions{
            1.0e-8, 1.0e-10, 30, false
        }
    };

    InitialState initial{
        {5.0e6, 5.0e6, 5.0e6},
        {10.0, 10.0},
        10.0
    };

    BoundarySchedule schedule{
        {5.0e6, 5.0e6},
        {10.0, 10.1}
    };

    const auto recorded =
        run_forward_scheduled_recorded(
            configuration, initial, schedule);

    require(recorded.summary.info == 0,
            "Fixture forward solve failed");

    auto reference = recorded_trajectory(recorded);

    FullSpaceProblemConfiguration cfg;
    cfg.reference = reference;
    cfg.parameters = configuration.parameters;

    const std::size_t width = 2*n+2;
    const std::size_t block = width+6;
    const std::size_t levels = reference.states.size();

    cfg.variables.assign(
        levels*block,
        ProblemVariable{false, 1.0});

    // Deliberately use distinct scales.
    cfg.variables[5] = {true, 1.0e6};
    cfg.variables[block+1] = {true, 10.0};
    cfg.variables[block+3] = {true, 2.0e6};
    cfg.variables[block+width+3] = {true, 3.0e6};

    const auto constraints =
        assemble_full_space_trajectory_constraints(
            reference, cfg.parameters);

    cfg.constraint_scales.assign(
        constraints.combined().size(), 1.0);

    cfg.measurements.resize(levels);

    for (std::size_t k = 0; k < levels; ++k) {
        cfg.measurements[k].observed =
            reference.reconstructed[k];

        cfg.measurements[k].available =
            {false, false, false, false};

        cfg.measurements[k].sigma =
            {1.0, 1.0, 1.0, 1.0};
    }

    const FullSpaceProblem problem(std::move(cfg));

    require(problem.free_variable_count() == 4,
            "Unexpected mapping fixture dimension");

    for (std::size_t j = 0; j < 4; ++j) {
        std::vector<double> direction(4, 0.0);
        direction[j] = 1.0;

        const auto tangent =
            problem.trajectory_direction(direction);

        double observed_scale = 0.0;

        if (j == 0)
            observed_scale = tangent.states[0][5];
        else if (j == 1)
            observed_scale = tangent.states[1][1];
        else if (j == 2)
            observed_scale = tangent.states[1][3];
        else
            observed_scale = tangent.reconstructed[1][1];

        const std::array<double, 4> expected{
            1.0e6, 10.0, 2.0e6, 3.0e6
        };

        near(
            observed_scale, expected[j],
            0.0, 0.0,
            "Incorrect physical coordinate scale");
    }

    std::cout
        << "PASS exact scaled-coordinate mapping\n";
}

void test_forecast_derivative()
{
    constexpr int n = 2;
    constexpr int steps = 10;

    const pipe_sim::TransientParameters parameters{
        100.0, 1.0, 0.65,
        0.5, 0.19634954084936207,
        288.15, 0.9, 518.3,
        1.1e-5, 1.0e-5
    };

    Configuration configuration{
        n,
        steps,
        parameters,
        pipe_sim::DemandRamp{
            5.0e6, 10.0, 10.0, 10.0
        },
        pipe_sim::NewtonOptions{
            1.0e-10, 1.0e-12, 40, false
        }
    };

    const InitialState initial{
        {5.0e6, 5.0e6, 5.0e6},
        {10.0, 10.0},
        10.0
    };

    BoundarySchedule schedule;
    schedule.inlet_pressure.assign(steps+1, 5.0e6);
    schedule.outlet_flow.assign(steps+1, 10.0);

    const auto nominal =
        run_forward_scheduled_recorded(
            configuration, initial, schedule);

    require(nominal.summary.info == 0,
            "Nominal forecast failed");

    const auto trajectory =
        recorded_trajectory(nominal);

    Eigen::VectorXd direction(2*n+2);
    direction <<
        0.0, 0.03, -0.02,
        1200.0, 0.04, -800.0;

    const auto tangent =
        propagate_forecast_tangent(
            trajectory, direction, parameters);

    require(
        tangent.states.size() ==
            static_cast<std::size_t>(steps+1),
        "Incomplete tangent trajectory");

    require(
        tangent.maximum_equation_residual < 1.0e-6,
        "Tangent linear solve residual too large");

    require(
        tangent.maximum_mass_balance_defect < 1.0e-8,
        "Linearized mass conservation failed");

    // Independent centered differences of the nonlinear
    // forecast, used only as a verification oracle.
    for (double h : {1.0e-2, 1.0e-3, 1.0e-4}) {
        auto plus = initial;
        auto minus = initial;

        const auto perturb = [&](InitialState& state,
                                 double sign) {
            state.pressure[0] += sign*h*direction[0];
            state.inlet_flow += sign*h*direction[1];

            for (int i = 0; i < n; ++i) {
                state.flow[i] +=
                    sign*h*direction[2*i+2];

                state.pressure[i+1] +=
                    sign*h*direction[2*i+3];
            }
        };

        perturb(plus, +1.0);
        perturb(minus, -1.0);

        const auto forward_plus =
            run_forward_scheduled(
                configuration, plus, schedule);

        const auto forward_minus =
            run_forward_scheduled(
                configuration, minus, schedule);

        require(
            forward_plus.info == 0 &&
            forward_minus.info == 0,
            "Perturbed nonlinear forecast failed");

        double maximum_pressure_error = 0.0;

        for (int k = 0; k <= steps; ++k) {
            const auto index =
                static_cast<std::size_t>(k);

            const double fd_pressure =
                (forward_plus.outlet_pressure[index] -
                 forward_minus.outlet_pressure[index]) /
                (2.0*h);

            const double exact_pressure =
                tangent.states[index][2*n+1];

            maximum_pressure_error =
                std::max(
                    maximum_pressure_error,
                    std::abs(fd_pressure-exact_pressure));

            near(
                fd_pressure,
                exact_pressure,
                2.0e-2,
                1.0e-3,
                "Outlet pressure directional derivative");

            const double fd_linepack =
                (forward_plus.linepack[index] -
                 forward_minus.linepack[index]) /
                (2.0*h);

            const auto& state = tangent.states[index];

            const double crho =
                pipe_sim::density_derivative_cz(
                    parameters.temperature,
                    parameters.z,
                    parameters.gas_constant);

            const double exact_linepack =
                parameters.area*parameters.dx*crho*
                (0.5*state[0] +
                 state[3] +
                 0.5*state[5]);

            near(
                fd_linepack,
                exact_linepack,
                1.0e-5,
                1.0e-3,
                "Linepack directional derivative");
        }

        std::cout
            << "h=" << h
            << " maximum pressure derivative error="
            << maximum_pressure_error
            << '\n';
    }

    std::cout
        << "maximum tangent equation residual="
        << tangent.maximum_equation_residual << '\n'
        << "maximum linearized mass defect="
        << tangent.maximum_mass_balance_defect << '\n'
        << "PASS independent tangent forecast validation\n";
}

} // namespace

int main()
{
    try {
        test_direction_mapping();
        test_forecast_derivative();
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E5.4 validation: "
            << error.what() << '\n';
        return 1;
    }
}
