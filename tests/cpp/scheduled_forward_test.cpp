#include "tpsf/scheduled_forward.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(double a, double b, double atol, double rtol,
          const std::string& message)
{
    require(std::isfinite(a) && std::isfinite(b),
            message + ": nonfinite value");

    require(std::abs(a - b) <=
                atol + rtol * std::max(std::abs(a), std::abs(b)),
            message + ": numerical disagreement");
}

void compare(const std::vector<double>& a,
             const std::vector<double>& b,
             double atol,
             double rtol,
             const std::string& name)
{
    require(a.size() == b.size(), name + ": size mismatch");

    for (std::size_t i = 0; i < a.size(); ++i)
        near(a[i], b[i], atol, rtol, name);
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

tpsf::Configuration configuration()
{
    constexpr int n = 10;
    constexpr double diameter = 0.5;
    constexpr double area = 0.19634954084936207;

    return {
        n, 40,
        pipe_sim::TransientParameters{
            100., 1., 0.65, diameter, area,
            288.15, 0.9, 518.3, 1.1e-5, 1.0e-5
        },
        pipe_sim::DemandRamp{
            5.0e6, 10.0, 12.0, 20.0
        },
        pipe_sim::NewtonOptions{
            1.0e-8, 1.0e-10, 30, false
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

tpsf::BoundarySchedule ramp_schedule(
    const tpsf::Configuration& cfg)
{
    tpsf::BoundarySchedule schedule;

    for (int k = 0; k <= cfg.steps; ++k) {
        const double time = k * cfg.parameters.dt;

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

double conservation_defect(
    const tpsf::Result& result,
    const tpsf::Configuration& cfg)
{
    double maximum = 0.0;

    const double dt = cfg.parameters.dt;
    const double theta = cfg.parameters.theta;

    for (int k = 1; k <= cfg.steps; ++k) {
        const double imbalance_new =
            result.inlet_flow[k] - result.outlet_flow[k];

        const double imbalance_old =
            result.inlet_flow[k-1] - result.outlet_flow[k-1];

        const double defect = std::abs(
            (result.linepack[k] - result.linepack[k-1]) / dt
            - (theta * imbalance_new +
               (1.0 - theta) * imbalance_old));

        maximum = std::max(maximum, defect);
    }

    return maximum;
}

void test_ramp_parity()
{
    const auto cfg = configuration();
    const auto initial = initial_state(cfg.cells);
    const auto schedule = ramp_schedule(cfg);

    const auto reference = tpsf::run_forward(cfg, initial);
    const auto scheduled =
        tpsf::run_forward_scheduled(cfg, initial, schedule);

    require(reference.info == 0, "Reference integration failed");
    require(scheduled.info == 0, "Scheduled integration failed");

    require(reference.completed_steps == cfg.steps,
            "Reference integration incomplete");

    require(scheduled.completed_steps == cfg.steps,
            "Scheduled integration incomplete");

    compare(reference.time, scheduled.time,
            0.0, 0.0, "Time");

    compare(reference.inlet_pressure,
            scheduled.inlet_pressure,
            0.0, 0.0, "Inlet pressure");

    compare(reference.inlet_flow,
            scheduled.inlet_flow,
            1.0e-10, 1.0e-12, "Inlet flow");

    compare(reference.outlet_flow,
            scheduled.outlet_flow,
            0.0, 0.0, "Outlet flow");

    compare(reference.outlet_pressure,
            scheduled.outlet_pressure,
            1.0e-7, 1.0e-12, "Outlet pressure");

    compare(reference.linepack,
            scheduled.linepack,
            1.0e-9, 1.0e-12, "Linepack");

    compare(reference.final_pressure,
            scheduled.final_pressure,
            1.0e-7, 1.0e-12, "Final pressure");

    compare(reference.final_flow,
            scheduled.final_flow,
            1.0e-10, 1.0e-12, "Final flow");

    near(reference.final_inlet_flow,
         scheduled.final_inlet_flow,
         1.0e-10, 1.0e-12, "Final inlet flow");

    std::cout << "PASS ramp parity\n";
}

void test_variable_boundaries()
{
    auto cfg = configuration();
    cfg.steps = 20;

    const auto initial = initial_state(cfg.cells);

    tpsf::BoundarySchedule schedule;

    for (int k = 0; k <= cfg.steps; ++k) {
        const double t = k * cfg.parameters.dt;

        schedule.inlet_pressure.push_back(
            5.0e6 + 1000.0 * std::sin(0.1 * t));

        schedule.outlet_flow.push_back(
            10.0 + 0.2 * std::sin(0.15 * t));
    }

    const auto result =
        tpsf::run_forward_scheduled(cfg, initial, schedule);

    require(result.info == 0, "Variable-boundary integration failed");

    require(result.completed_steps == cfg.steps,
            "Variable-boundary integration incomplete");

    compare(result.inlet_pressure,
            schedule.inlet_pressure,
            0.0, 0.0, "Scheduled inlet pressure");

    compare(result.outlet_flow,
            schedule.outlet_flow,
            0.0, 0.0, "Scheduled outlet flow");

    const double defect = conservation_defect(result, cfg);

    near(defect, result.max_mass_defect,
         1.0e-12, 1.0e-12,
         "Reported conservation defect");

    require(defect < 1.0e-5,
            "Mass conservation defect exceeds acceptance threshold");

    std::cout
        << "PASS variable boundaries; max mass defect="
        << defect << " kg/s\n";
}

void test_invalid_schedules()
{
    const auto cfg = configuration();
    const auto initial = initial_state(cfg.cells);
    const auto valid = ramp_schedule(cfg);

    auto invalid = valid;
    invalid.outlet_flow.pop_back();

    expect_invalid([&] {
        tpsf::run_forward_scheduled(cfg, initial, invalid);
    });

    invalid = valid;
    invalid.inlet_pressure[1] =
        std::numeric_limits<double>::quiet_NaN();

    expect_invalid([&] {
        tpsf::run_forward_scheduled(cfg, initial, invalid);
    });

    invalid = valid;
    invalid.inlet_pressure[0] += 1000.0;

    expect_invalid([&] {
        tpsf::run_forward_scheduled(cfg, initial, invalid);
    });

    std::cout << "PASS invalid schedule rejection\n";
}

} // namespace

int main()
{
    try {
        test_ramp_parity();
        test_variable_boundaries();
        test_invalid_schedules();

        std::cout << "PASS E2.3 scheduled integration\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FAIL E2.3: " << e.what() << '\n';
        return 1;
    }
}
