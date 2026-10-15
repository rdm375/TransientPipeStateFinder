#include "tpsf/observations.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

template<typename Function>
void expect_invalid(Function function)
{
    bool rejected = false;

    try {
        function();
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    require(rejected, "Invalid observation record was accepted");
}

} // namespace

int main()
{
    try {
        using namespace tpsf;

        BoundaryObservationRecord record;

        record.measured.inlet_pressure = {
            8.0e6, 8.0e6, 8.0e6
        };

        record.measured.inlet_flow = {
            100.0, 101.0, 102.0
        };

        record.measured.outlet_pressure = {
            7.9e6, 7.89e6, 7.88e6
        };

        record.measured.outlet_flow = {
            100.0, 105.0, 105.0
        };

        record.scales.inlet_flow = 2.0;
        record.scales.outlet_pressure = 5000.0;

        record.include_initial = false;

        validate_boundary_observations(record, 3);
        validate_current_boundary_roles(record.roles);

        const auto observations =
            make_current_observations(record, 3);

        require(observations.inlet_flow.size() == 3,
                "Incorrect inlet-flow observation count");

        require(observations.outlet_pressure.size() == 3,
                "Incorrect outlet-pressure observation count");

        require(observations.inlet_flow[1] == 101.0,
                "Incorrect inlet-flow observation");

        require(observations.outlet_pressure[2] == 7.88e6,
                "Incorrect outlet-pressure observation");

        require(observations.sigma_flow == 2.0,
                "Incorrect flow observation scale");

        require(observations.sigma_pressure == 5000.0,
                "Incorrect pressure observation scale");

        require(!observations.include_initial,
                "Incorrect initial-observation setting");

        // All four channels must have consistent dimensions.
        auto invalid = record;
        invalid.measured.outlet_flow.pop_back();

        expect_invalid([&] {
            validate_boundary_observations(invalid, 3);
        });

        // Measurement scales must be positive.
        invalid = record;
        invalid.scales.inlet_flow = 0.0;

        expect_invalid([&] {
            validate_boundary_observations(invalid, 3);
        });

        // Nonfinite observations are forbidden.
        invalid = record;
        invalid.measured.outlet_pressure[1] =
            std::numeric_limits<double>::quiet_NaN();

        expect_invalid([&] {
            validate_boundary_observations(invalid, 3);
        });

        // The current formulation cannot use prescribed
        // inlet pressure as an independent residual.
        invalid = record;
        invalid.roles.inlet_pressure =
            MeasurementRole::ObjectiveResidual;

        expect_invalid([&] {
            make_current_observations(invalid, 3);
        });

        std::cout << "PASS boundary observation model\n";

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
