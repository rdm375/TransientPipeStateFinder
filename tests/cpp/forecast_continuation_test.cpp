#include "tpsf/forecast_continuation.hpp"

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
    double a,
    double b,
    double atol,
    double rtol,
    const std::string& message)
{
    require(
        std::isfinite(a) && std::isfinite(b),
        message + ": nonfinite value");

    require(
        std::abs(a-b) <=
            atol + rtol*std::max(std::abs(a), std::abs(b)),
        message + ": numerical mismatch");
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

    require(rejected, "Invalid continuation was accepted");
}

tpsf::FullSpaceTrajectory reconstructed_from_recording(
    const tpsf::RecordedForwardResult& recorded)
{
    tpsf::FullSpaceTrajectory trajectory;

    trajectory.n = recorded.trajectory.cells;

    for (const auto& snapshot :
         recorded.trajectory.snapshots) {
        trajectory.states.push_back(snapshot.state);
        trajectory.inputs.push_back(snapshot.inputs);
        trajectory.reconstructed.push_back(
            snapshot.physical_channels);
    }

    return trajectory;
}

void test_split_continuation()
{
    constexpr int n = 10;
    constexpr int total_steps = 20;
    constexpr int split = 8;

    const tpsf::Configuration configuration{
        n,
        total_steps,
        pipe_sim::TransientParameters{
            100., 1., 0.65,
            0.5, 0.19634954084936207,
            288.15, 0.9, 518.3,
            1.1e-5, 1.0e-5
        },
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

    for (int k = 0; k <= total_steps; ++k) {
        const double t = k * configuration.parameters.dt;

        schedule.inlet_pressure.push_back(
            5.0e6 + 1000.0*std::sin(0.1*t));

        schedule.outlet_flow.push_back(
            10.0 + 0.2*std::sin(0.15*t));
    }

    // Initial inlet pressure must match initial state.
    schedule.inlet_pressure[0] = initial.pressure[0];

    const auto uninterrupted =
        tpsf::run_forward_scheduled_recorded(
            configuration, initial, schedule);

    require(
        uninterrupted.summary.info == 0 &&
        uninterrupted.summary.completed_steps == total_steps,
        "Uninterrupted integration failed");

    auto reconstruction_configuration = configuration;
    reconstruction_configuration.steps = split;

    tpsf::BoundarySchedule past;

    past.inlet_pressure.assign(
        schedule.inlet_pressure.begin(),
        schedule.inlet_pressure.begin()+split+1);

    past.outlet_flow.assign(
        schedule.outlet_flow.begin(),
        schedule.outlet_flow.begin()+split+1);

    const auto past_recording =
        tpsf::run_forward_scheduled_recorded(
            reconstruction_configuration,
            initial,
            past);

    require(
        past_recording.summary.info == 0 &&
        past_recording.summary.completed_steps == split,
        "Past integration failed");

    const auto reconstructed =
        reconstructed_from_recording(past_recording);

    auto forecast_configuration = configuration;
    forecast_configuration.steps = total_steps-split;

    tpsf::BoundarySchedule future;

    future.inlet_pressure.assign(
        schedule.inlet_pressure.begin()+split,
        schedule.inlet_pressure.end());

    future.outlet_flow.assign(
        schedule.outlet_flow.begin()+split,
        schedule.outlet_flow.end());

    const auto continued =
        tpsf::forecast_from_full_space(
            forecast_configuration,
            reconstructed,
            future);

    require(
        continued.summary.info == 0 &&
        continued.summary.completed_steps ==
            forecast_configuration.steps,
        "Forecast continuation failed");

    for (int k = 0;
         k <= forecast_configuration.steps;
         ++k) {
        const auto original =
            static_cast<std::size_t>(split+k);

        const auto local =
            static_cast<std::size_t>(k);

        near(
            continued.summary.outlet_pressure[local],
            uninterrupted.summary.outlet_pressure[original],
            1.0e-7, 1.0e-12,
            "Outlet pressure");

        near(
            continued.summary.inlet_flow[local],
            uninterrupted.summary.inlet_flow[original],
            1.0e-10, 1.0e-12,
            "Inlet flow");

        near(
            continued.summary.linepack[local],
            uninterrupted.summary.linepack[original],
            1.0e-9, 1.0e-12,
            "Linepack");
    }

    require(
        continued.summary.max_mass_defect < 1.0e-5,
        "Forecast conservation defect too large");

    auto invalid = future;
    invalid.outlet_flow[0] += 1.0;

    expect_invalid([&] {
        (void)tpsf::forecast_from_full_space(
            forecast_configuration,
            reconstructed,
            invalid);
    });

    invalid = future;
    invalid.inlet_pressure[0] += 1000.0;

    expect_invalid([&] {
        (void)tpsf::forecast_from_full_space(
            forecast_configuration,
            reconstructed,
            invalid);
    });

    std::cout
        << "PASS split forecast continuation\n"
        << "max mass defect: "
        << continued.summary.max_mass_defect
        << " kg/s\n";
}

} // namespace

int main()
{
    try {
        test_split_continuation();
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr
            << "FAIL forecast continuation: "
            << e.what() << '\n';
        return 1;
    }
}
