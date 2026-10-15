#include "tpsf/objective.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected,
          double atol, double rtol,
          const char* message)
{
    require(std::isfinite(actual) && std::isfinite(expected),
            "Nonfinite comparison");

    if (std::abs(actual - expected) >
        atol + rtol * std::max(std::abs(actual), std::abs(expected)))
        throw std::runtime_error(message);
}

template<class Function>
void expect_invalid(Function function)
{
    bool rejected = false;

    try {
        function();
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    require(rejected, "Invalid input was accepted");
}

tpsf::ObjectiveProblem make_problem()
{
    constexpr int n = 10;
    constexpr int steps = 20;
    constexpr double diameter = 0.5;
    constexpr double area = 0.19634954084936207;

    tpsf::ObjectiveProblem problem;

    problem.configuration = {
        n, steps,
        pipe_sim::TransientParameters{
            100., 1., 0.65, diameter, area,
            288.15, 0.9, 518.3, 1.1e-5, 1.0e-5
        },
        pipe_sim::DemandRamp{
            5.0e6, 10.0, 10.0, 0.0
        },
        pipe_sim::NewtonOptions{
            1.0e-8, 1.0e-10, 30, false
        }
    };

    problem.baseline_initial = {
        std::vector<double>(n+1, 5.0e6),
        std::vector<double>(n, 10.0),
        10.0
    };

    auto& obs = problem.observations;

    obs.include_initial = true;

    for (int k = 0; k <= steps; ++k) {
        obs.time.push_back(static_cast<double>(k));

        obs.measured.inlet_pressure.push_back(5.0e6);
        obs.measured.inlet_flow.push_back(10.0);
        obs.measured.outlet_pressure.push_back(5.0e6);
        obs.measured.outlet_flow.push_back(10.0);
    }

    obs.scales.inlet_pressure = 1000.0;
    obs.scales.inlet_flow = 1.0;
    obs.scales.outlet_pressure = 1000.0;
    obs.scales.outlet_flow = 1.0;

    return problem;
}

void test_objective()
{
    auto problem = make_problem();

    const auto layout =
        tpsf::objective_layout(problem.configuration);

    require(layout.initial_size == 21,
            "Incorrect initial-state dimension");

    require(layout.boundary_size == 21,
            "Incorrect boundary dimension");

    require(layout.total_size == 63,
            "Incorrect optimization dimension");

    auto parameters =
        tpsf::zero_objective_parameters(problem);

    const auto baseline =
        tpsf::evaluate_objective(problem, parameters);

    require(baseline.forward.info == 0,
            "Baseline integration failed");

    require(std::isfinite(baseline.total_cost),
            "Nonfinite baseline cost");

    // Inlet-pressure correction must enter the cost.
    parameters[layout.inlet_pressure_offset + 10] = 0.1;

    const auto pressure_case =
        tpsf::evaluate_objective(problem, parameters);

    near(pressure_case.residuals.inlet_pressure[10],
         0.1, 1.0e-12, 1.0e-12,
         "Incorrect inlet-pressure residual");

    near(pressure_case.inlet_pressure_cost,
         0.005, 1.0e-12, 1.0e-12,
         "Incorrect inlet-pressure cost");

    require(
        pressure_case.predicted.inlet_pressure[10] !=
        baseline.predicted.inlet_pressure[10],
        "Inlet-pressure correction did not change prediction");

    // Outlet-flow correction must enter the cost.
    parameters =
        tpsf::zero_objective_parameters(problem);

    parameters[layout.outlet_flow_offset + 10] = 0.1;

    const auto flow_case =
        tpsf::evaluate_objective(problem, parameters);

    near(flow_case.residuals.outlet_flow[10],
         0.1, 1.0e-12, 1.0e-12,
         "Incorrect outlet-flow residual");

    near(flow_case.outlet_flow_cost,
         0.005, 1.0e-12, 1.0e-12,
         "Incorrect outlet-flow cost");

    // Verify initial inlet-pressure consistency when alpha_0
    // is perturbed. This must not cause a schedule mismatch.
    parameters =
        tpsf::zero_objective_parameters(problem);

    parameters[layout.inlet_pressure_offset] = 0.1;

    const auto initial_boundary_case =
        tpsf::evaluate_objective(problem, parameters);

    near(initial_boundary_case.residuals.inlet_pressure[0],
         0.1, 1.0e-12, 1.0e-12,
         "Incorrect initial inlet-pressure residual");

    // All four cost contributions must sum to the total.
    const double sum =
        flow_case.inlet_pressure_cost +
        flow_case.inlet_flow_cost +
        flow_case.outlet_pressure_cost +
        flow_case.outlet_flow_cost;

    near(sum, flow_case.measurement_cost,
         1.0e-12, 1.0e-12,
         "Incorrect measurement-cost decomposition");

    near(flow_case.total_cost,
         flow_case.measurement_cost,
         1.0e-12, 1.0e-12,
         "Unexpected regularization cost");

    // Invalid optimization-vector dimension.
    parameters.pop_back();

    expect_invalid([&] {
        tpsf::evaluate_objective(problem, parameters);
    });

    // Nonfinite optimization parameter.
    parameters =
        tpsf::zero_objective_parameters(problem);

    parameters[0] =
        std::numeric_limits<double>::quiet_NaN();

    expect_invalid([&] {
        tpsf::evaluate_objective(problem, parameters);
    });

    // Observation timestamps must match the integration grid.
    auto invalid_problem = make_problem();
    invalid_problem.observations.time[10] += 0.1;

    expect_invalid([&] {
        tpsf::zero_objective_parameters(invalid_problem);
    });

    std::cout
        << "PASS four-channel objective evaluation\n";
}

} // namespace

int main()
{
    try {
        test_objective();
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FAIL E2.4a: " << e.what() << '\n';
        return 1;
    }
}
