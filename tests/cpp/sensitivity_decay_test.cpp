#include "tpsf/sensitivity_decay.hpp"

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

void test_monotonic_decay()
{
    const std::vector<double> time{0, 1, 2, 3, 4};
    const std::vector<double> u{10, 5, 2, 0.5, 0.1};

    SensitivityDecayOptions options;
    options.absolute_threshold = 1.0;
    options.relative_threshold = 0.1;

    const auto result =
        assess_sensitivity_decay(time, u, options);

    require(!result.nonmonotonic,
            "Monotonic decay incorrectly classified");

    require(!result.delayed_peak,
            "Unexpected delayed peak");

    require(result.absolute.time.has_value() &&
            *result.absolute.time == 3.0,
            "Incorrect absolute attenuation time");

    require(result.relative.time.has_value() &&
            *result.relative.time == 3.0,
            "Incorrect relative attenuation time");
}

void test_rebound()
{
    const std::vector<double> time{0, 1, 2, 3, 4};
    const std::vector<double> u{10, 0.1, 5, 0.2, 0.1};

    SensitivityDecayOptions options;
    options.absolute_threshold = 1.0;

    const auto result =
        assess_sensitivity_decay(time, u, options);

    require(result.nonmonotonic,
            "Failed to detect sensitivity rebound");

    require(result.absolute.time.has_value() &&
            *result.absolute.time == 3.0,
            "Temporary crossing incorrectly accepted");
}

void test_delayed_peak()
{
    const std::vector<double> time{0, 1, 2, 3};
    const std::vector<double> u{1, 3, 2, 0.1};

    SensitivityDecayOptions options;
    options.absolute_threshold = 0.5;

    const auto result =
        assess_sensitivity_decay(time, u, options);

    require(result.delayed_peak,
            "Failed to detect delayed peak");

    require(result.peak_time == 1.0,
            "Incorrect delayed peak time");

    require(result.absolute.time.has_value() &&
            *result.absolute.time == 3.0,
            "Incorrect attenuation after delayed peak");
}

void test_not_attained()
{
    const std::vector<double> time{0, 1, 2};
    const std::vector<double> u{10, 5, 2};

    SensitivityDecayOptions options;
    options.absolute_threshold = 1.0;

    const auto result =
        assess_sensitivity_decay(time, u, options);

    require(result.absolute.status ==
                AttenuationStatus::not_attained,
            "Incorrectly reported attenuation");

    require(!result.absolute.time.has_value(),
            "Unexpected attenuation time");
}

void test_zero_reference()
{
    const std::vector<double> time{0, 1, 2};
    const std::vector<double> u{0, 2, 0};

    const auto result =
        assess_sensitivity_decay(time, u);

    require(result.relative.status ==
                AttenuationStatus::undefined_reference,
            "Zero initial reference must be undefined");

    require(result.delayed_peak,
            "Expected delayed peak");
}

void test_forecast_adapter()
{
    ForecastDirectionalUncertainty forecast;

    forecast.time = {0, 1, 2};
    forecast.measurement_unit_radius = {3, 2, 0};
    forecast.unresolved_unit_radius = {0, 1, 0};

    SensitivityDecayOptions options;
    options.absolute_threshold = 0.5;

    const auto result =
        assess_forecast_sensitivity_decay(forecast, options);

    require(result.measurement.absolute.time.has_value() &&
            *result.measurement.absolute.time == 2,
            "Incorrect measurement attenuation");

    require(result.unresolved.relative.status ==
                AttenuationStatus::undefined_reference,
            "Expected undefined unresolved relative reference");
}


void test_peak_relative_attenuation()
{
    SensitivityDecayOptions options;
    options.relative_threshold = 0.1;

    // Monotonic decay: initial and peak references coincide.
    {
        const std::vector<double> time{0, 1, 2, 3, 4};
        const std::vector<double> u{10, 5, 2, 0.5, 0.1};

        const auto result =
            assess_sensitivity_decay(time, u, options);

        require(result.peak_relative.time.has_value() &&
                *result.peak_relative.time == 3.0,
                "Incorrect monotonic peak-relative time");

        require(result.relative.time ==
                    result.peak_relative.time,
                "Initial and peak references should coincide");
    }

    // Initially zero sensitivity with a delayed peak.
    {
        const std::vector<double> time{0, 1, 2, 3, 4};
        const std::vector<double> u{0, 10, 5, 1, 0.1};

        const auto result =
            assess_sensitivity_decay(time, u, options);

        require(result.relative.status ==
                    AttenuationStatus::undefined_reference,
                "Zero initial reference must remain undefined");

        require(result.peak_relative.status ==
                    AttenuationStatus::attained,
                "Expected peak-relative attenuation");

        require(result.peak_relative.time.has_value() &&
                *result.peak_relative.time == 3.0,
                "Incorrect delayed peak-relative time");
    }

    // A temporary crossing followed by a rebound must not
    // be accepted as persistent attenuation.
    {
        const std::vector<double> time{0, 1, 2, 3, 4};
        const std::vector<double> u{10, 0.1, 5, 0.2, 0.1};

        const auto result =
            assess_sensitivity_decay(time, u, options);

        require(result.peak_relative.time.has_value() &&
                *result.peak_relative.time == 3.0,
                "Peak-relative attenuation ignored rebound");
    }

    // An identically zero history has no positive reference
    // magnitude from which to define a relative threshold.
    {
        const std::vector<double> time{0, 1, 2};
        const std::vector<double> u{0, 0, 0};

        const auto result =
            assess_sensitivity_decay(time, u, options);

        require(result.peak_relative.status ==
                    AttenuationStatus::undefined_reference,
                "Zero peak reference must be undefined");

        require(result.absolute.status ==
                    AttenuationStatus::attained,
                "Zero sensitivity must satisfy absolute threshold");
    }

    std::cout
        << "PASS peak-relative attenuation regression\n";
}

} // namespace

int main()
{
    try {
        test_monotonic_decay();
        test_rebound();
        test_delayed_peak();
        test_not_attained();
        test_zero_reference();
        test_forecast_adapter();
        test_peak_relative_attenuation();

        std::cout
            << "PASS E5.6a sensitivity decay diagnostics\n";

        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL E5.6a: "
                  << error.what() << '\n';
        return 1;
    }
}
