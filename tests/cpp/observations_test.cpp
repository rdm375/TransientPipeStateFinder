#include "tpsf/boundary_residuals.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message)
{
    if (!condition)
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

    require(rejected, "Expected invalid_argument");
}

} // namespace

int main()
{
    try {
        using namespace tpsf;

        BoundaryObservationRecord observations;

        observations.time = {0.0, 10.0, 20.0};

        observations.measured.inlet_pressure =
            {8.0e6, 8.0e6, 8.0e6};

        observations.measured.inlet_flow =
            {100.0, 101.0, 102.0};

        observations.measured.outlet_pressure =
            {7.9e6, 7.89e6, 7.88e6};

        observations.measured.outlet_flow =
            {100.0, 105.0, 105.0};

        observations.scales.inlet_pressure = 1000.0;
        observations.scales.inlet_flow = 2.0;
        observations.scales.outlet_pressure = 1000.0;
        observations.scales.outlet_flow = 2.0;

        validate_boundary_observations(observations);

        BoundaryMeasurements predicted = observations.measured;

        // One normalized residual of +1 per channel.
        predicted.inlet_pressure[1] += 1000.0;
        predicted.inlet_flow[1] += 2.0;
        predicted.outlet_pressure[1] += 1000.0;
        predicted.outlet_flow[1] += 2.0;

        const auto residuals =
            compute_boundary_residuals(predicted, observations);

        const double cost =
            boundary_least_squares_cost(residuals);

        require(std::abs(cost - 2.0) < 1.0e-12,
                "Incorrect four-channel cost");

        require(residuals.inlet_pressure[1] == 1.0,
                "Incorrect inlet-pressure residual");

        require(residuals.inlet_flow[1] == 1.0,
                "Incorrect inlet-flow residual");

        require(residuals.outlet_pressure[1] == 1.0,
                "Incorrect outlet-pressure residual");

        require(residuals.outlet_flow[1] == 1.0,
                "Incorrect outlet-flow residual");

        // Zero residual must produce zero cost.
        const auto zero =
            compute_boundary_residuals(
                observations.measured, observations);

        require(boundary_least_squares_cost(zero) == 0.0,
                "Identical data must have zero cost");

        // Reject non-monotonic timestamps.
        auto invalid = observations;
        invalid.time[2] = invalid.time[1];

        expect_invalid([&] {
            validate_boundary_observations(invalid);
        });

        // Reject invalid uncertainties.
        invalid = observations;
        invalid.scales.inlet_pressure = 0.0;

        expect_invalid([&] {
            validate_boundary_observations(invalid);
        });

        std::cout << "PASS four-channel observation model\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
