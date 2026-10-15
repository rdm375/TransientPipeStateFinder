#include "tpsf/full_space_measurement_objective.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
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

void near(double actual, double expected,
          double atol, double rtol,
          const std::string& message)
{
    require(std::isfinite(actual) &&
            std::isfinite(expected),
            message + ": nonfinite");

    require(
        std::abs(actual-expected) <=
        atol + rtol*std::max(
            std::abs(actual), std::abs(expected)),
        message + ": mismatch");
}

template<class F>
void expect_invalid(F&& f)
{
    bool caught = false;

    try {
        f();
    }
    catch (const std::invalid_argument&) {
        caught = true;
    }

    require(caught, "Expected invalid_argument");
}

struct Fixture {
    std::vector<BoundaryVector> reconstructed;
    std::vector<BoundaryMeasurement> measurements;

    Fixture()
    {
        reconstructed = {
            BoundaryVector{5.0e6, 10.0, 4.9e6, 11.0},
            BoundaryVector{5.01e6, 10.5, 4.91e6, 11.5},
            BoundaryVector{5.02e6, 11.0, 4.92e6, 12.0}
        };

        measurements.resize(3);

        for (std::size_t k = 0; k < 3; ++k) {
            for (std::size_t c = 0; c < 4; ++c) {
                const double sigma =
                    (c == 0 || c == 2) ? 100.0 : 0.2;

                measurements[k].sigma[c] = sigma;

                // Exactly one normalized residual unit
                // on each available channel.
                measurements[k].observed[c] =
                    reconstructed[k][c] - sigma;

                measurements[k].available[c] = true;
            }
        }
    }
};

void test_all_channels()
{
    Fixture f;

    const auto result =
        evaluate_full_space_measurement_objective(
            f.reconstructed, f.measurements);

    require(result.measurement_count == 12,
            "Incorrect measurement count");

    near(result.value, 6.0, 1e-10, 1e-12,
         "All-channel objective");

    require(result.normalized_rms.has_value(),
            "Missing normalized RMS");

    require(result.maximum_normalized_error.has_value(),
            "Missing maximum error");

    near(*result.normalized_rms, 1.0,
         1e-10, 1e-12, "Normalized RMS");

    near(*result.maximum_normalized_error, 1.0,
         1e-10, 1e-12, "Maximum normalized error");

    for (std::size_t k = 0; k < 3; ++k) {
        for (std::size_t c = 0; c < 4; ++c) {
            near(result.gradient[k][c],
                 1.0/f.measurements[k].sigma[c],
                 1e-10, 1e-10,
                 "Analytic gradient");
        }
    }
}

void test_missing_channels()
{
    Fixture f;

    f.measurements[0].available[0] = false;
    f.measurements[1].available[1] = false;
    f.measurements[2].available[2] = false;

    const double nan =
        std::numeric_limits<double>::quiet_NaN();

    f.measurements[0].observed[0] = nan;
    f.measurements[0].sigma[0] = nan;

    f.measurements[1].observed[1] = nan;
    f.measurements[1].sigma[1] = 0.0;

    f.measurements[2].observed[2] = nan;
    f.measurements[2].sigma[2] = -1.0;

    const auto result =
        evaluate_full_space_measurement_objective(
            f.reconstructed, f.measurements);

    require(result.measurement_count == 9,
            "Incorrect masked measurement count");

    near(result.value, 4.5, 1e-10, 1e-12,
         "Masked objective");

    require(result.gradient[0][0] == 0.0,
            "Masked pressure gradient nonzero");

    require(result.gradient[1][1] == 0.0,
            "Masked flow gradient nonzero");

    require(result.gradient[2][2] == 0.0,
            "Masked outlet gradient nonzero");
}

void test_no_measurements()
{
    Fixture f;

    for (auto& measurement : f.measurements)
        measurement.available.fill(false);

    const auto result =
        evaluate_full_space_measurement_objective(
            f.reconstructed, f.measurements);

    require(result.value == 0.0,
            "Empty objective is nonzero");

    require(result.measurement_count == 0,
            "Empty measurement count nonzero");

    require(!result.normalized_rms.has_value(),
            "Undefined RMS reported as valid");

    require(!result.maximum_normalized_error.has_value(),
            "Undefined maximum reported as valid");

    for (const auto& gradient : result.gradient)
        for (double component : gradient)
            require(component == 0.0,
                    "Empty measurement gradient nonzero");
}

void test_gradient()
{
    Fixture f;

    f.measurements[0].available[1] = false;
    f.measurements[2].available[3] = false;

    const auto result =
        evaluate_full_space_measurement_objective(
            f.reconstructed, f.measurements);

    // Different finite-difference steps for pressure
    // and mass-flow channels.
    for (std::size_t k = 0; k < f.reconstructed.size(); ++k) {
        for (std::size_t c = 0; c < 4; ++c) {
            const double h =
                (c == 0 || c == 2) ? 0.01 : 1e-5;

            auto plus = f.reconstructed;
            auto minus = f.reconstructed;

            plus[k][c] += h;
            minus[k][c] -= h;

            const double jp =
                evaluate_full_space_measurement_objective(
                    plus, f.measurements).value;

            const double jm =
                evaluate_full_space_measurement_objective(
                    minus, f.measurements).value;

            const double numerical = (jp-jm)/(2*h);

            near(numerical, result.gradient[k][c],
                 1e-7, 2e-7,
                 "Finite-difference gradient");
        }
    }
}

void test_validation()
{
    Fixture f;

    auto shorter = f.measurements;
    shorter.pop_back();

    expect_invalid([&] {
        evaluate_full_space_measurement_objective(
            f.reconstructed, shorter);
    });

    f.measurements[0].sigma[0] = 0.0;

    expect_invalid([&] {
        evaluate_full_space_measurement_objective(
            f.reconstructed, f.measurements);
    });

    f.measurements[0].sigma[0] = -1.0;

    expect_invalid([&] {
        evaluate_full_space_measurement_objective(
            f.reconstructed, f.measurements);
    });

    f.measurements[0].sigma[0] =
        std::numeric_limits<double>::infinity();

    expect_invalid([&] {
        evaluate_full_space_measurement_objective(
            f.reconstructed, f.measurements);
    });

    f.measurements[0].sigma[0] = 100.0;
    f.measurements[0].observed[0] =
        std::numeric_limits<double>::quiet_NaN();

    expect_invalid([&] {
        evaluate_full_space_measurement_objective(
            f.reconstructed, f.measurements);
    });

    f.measurements[0].observed[0] = 5.0e6;
    f.reconstructed[0][0] =
        std::numeric_limits<double>::infinity();

    expect_invalid([&] {
        evaluate_full_space_measurement_objective(
            f.reconstructed, f.measurements);
    });
}

void test_exact_zero()
{
    Fixture f;

    for (std::size_t k = 0; k < f.measurements.size(); ++k)
        for (std::size_t c = 0; c < 4; ++c)
            f.measurements[k].observed[c] =
                f.reconstructed[k][c];

    const auto result =
        evaluate_full_space_measurement_objective(
            f.reconstructed, f.measurements);

    require(result.value == 0.0,
            "Exact measurements produce nonzero objective");

    for (const auto& gradient : result.gradient)
        for (double component : gradient)
            require(component == 0.0,
                    "Exact measurements produce nonzero gradient");
}

} // namespace

int main()
{
    try {
        test_all_channels();
        test_missing_channels();
        test_no_measurements();
        test_gradient();
        test_validation();
        test_exact_zero();

        std::cout
            << "PASS E2.4b.9 measurement objective\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E2.4b.9: "
            << error.what() << '\n';

        return 1;
    }
}
