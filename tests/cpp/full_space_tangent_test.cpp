#include "tpsf/full_space_tangent.hpp"

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

struct Fixture {
    static constexpr int n = 5;

    pipe_sim::TransientParameters parameters{
        100.0, 1.0, 0.65,
        0.5, 0.19634954084936207,
        288.15, 0.9, 518.3,
        1.1e-5, 1.0e-5
    };

    std::vector<double> old_state{
        5.00e6, 10.2, 10.1, 4.98e6,
        10.4, 4.96e6, 10.7, 4.94e6,
        11.0, 4.92e6, 11.3, 4.90e6
    };

    std::vector<double> new_state{
        5.01e6, 10.5, 10.3, 4.99e6,
        10.6, 4.97e6, 10.9, 4.95e6,
        11.2, 4.93e6, 11.5, 4.91e6
    };

    tpsf::PhysicalBoundaryInputs old_inputs{
        5.00e6, 11.4
    };

    tpsf::PhysicalBoundaryInputs new_inputs{
        5.01e6, 11.6
    };

    tpsf::BoundaryVector reconstructed{
        5.01e6, 10.5, 4.91e6, 11.6
    };
};

struct Direction {
    std::vector<double> old_state;
    std::vector<double> new_state;

    tpsf::PhysicalBoundaryInputs old_inputs{};
    tpsf::PhysicalBoundaryInputs new_inputs{};

    tpsf::BoundaryVector reconstructed{};
};

Direction make_direction(int group)
{
    Direction d;

    d.old_state.assign(2*Fixture::n+2, 0.0);
    d.new_state.assign(2*Fixture::n+2, 0.0);

    // Pressure perturbations are measured in Pa.
    // Flow perturbations are measured in kg/s.
    //
    // These directions are scaled so that h around
    // 1e-3 to 1e-2 produces resolvable perturbations.

    if (group == 0 || group == 5) {
        for (std::size_t i = 0; i < d.old_state.size(); ++i) {
            d.old_state[i] =
                (i == 0 || (i >= 3 && i % 2 == 1))
                    ? 1.0e5*(1.0 + 0.1*i)
                    : 2.0*(1.0 + 0.05*i);
        }
    }

    if (group == 1 || group == 5) {
        for (std::size_t i = 0; i < d.new_state.size(); ++i) {
            d.new_state[i] =
                (i == 0 || (i >= 3 && i % 2 == 1))
                    ? -8.0e4*(1.0 + 0.08*i)
                    : 1.5*(1.0 + 0.04*i);
        }
    }

    if (group == 2 || group == 5) {
        d.old_inputs = {1.2e5, -2.0};
    }

    if (group == 3 || group == 5) {
        d.new_inputs = {-1.1e5, 1.7};
    }

    if (group == 4 || group == 5) {
        d.reconstructed = {
            9.0e4, -1.3, -7.0e4, 2.1
        };
    }

    return d;
}

std::vector<double> evaluate(
    const Fixture& f,
    const Direction& d,
    double h)
{
    auto old_state = f.old_state;
    auto new_state = f.new_state;

    for (std::size_t i = 0; i < old_state.size(); ++i)
        old_state[i] += h*d.old_state[i];

    for (std::size_t i = 0; i < new_state.size(); ++i)
        new_state[i] += h*d.new_state[i];

    auto old_inputs = f.old_inputs;
    auto new_inputs = f.new_inputs;
    auto reconstructed = f.reconstructed;

    old_inputs.inlet_pressure +=
        h*d.old_inputs.inlet_pressure;

    old_inputs.outlet_flow +=
        h*d.old_inputs.outlet_flow;

    new_inputs.inlet_pressure +=
        h*d.new_inputs.inlet_pressure;

    new_inputs.outlet_flow +=
        h*d.new_inputs.outlet_flow;

    for (std::size_t i = 0; i < reconstructed.size(); ++i)
        reconstructed[i] += h*d.reconstructed[i];

    return tpsf::assemble_full_space_step_constraints(
        Fixture::n,
        old_state,
        new_state,
        old_inputs,
        new_inputs,
        reconstructed,
        f.parameters
    ).combined();
}

void test_direction(int group)
{
    const Fixture f;
    const auto d = make_direction(group);

    const auto analytic =
        tpsf::apply_full_space_step_jvp(
            Fixture::n,
            f.old_state,
            f.new_state,
            d.old_state,
            d.new_state,
            d.old_inputs,
            d.new_inputs,
            d.reconstructed,
            f.parameters
        ).combined();

    const std::vector<double> steps{
        1.0e-2, 5.0e-3, 2.5e-3, 1.25e-3
    };

    std::vector<double> errors;

    for (double h : steps) {
        const auto plus = evaluate(f, d, h);
        const auto minus = evaluate(f, d, -h);

        require(plus.size() == analytic.size(),
                "Constraint size mismatch");

        double max_error = 0.0;

        for (std::size_t i = 0; i < analytic.size(); ++i) {
            const double numerical =
                (plus[i]-minus[i])/(2.0*h);

            const double scale =
                std::max(1.0, std::abs(analytic[i]));

            const double error =
                std::abs(numerical-analytic[i])/scale;

            require(std::isfinite(error),
                    "Nonfinite derivative error");

            max_error = std::max(max_error, error);
        }

        errors.push_back(max_error);
    }

    std::cout << "Direction group " << group << ": ";

    for (double error : errors)
        std::cout << error << " ";

    std::cout << '\n';

    // All derivative groups must meet the final accuracy
    // criterion at the smallest tested perturbation.
    require(
        errors.back() < 1.0e-8,
        "Tangent-linear derivative accuracy failed");

    // The nonlinear state derivatives should exhibit
    // second-order central-difference convergence.
    //
    // A factor of four is expected when h is halved.
    // Allow modest numerical variation.
    if (group == 0 || group == 1 || group == 5) {
        require(
            errors[0] > 3.0*errors[1] &&
            errors[1] > 3.0*errors[2],
            "Second-order derivative convergence failed");
    }
}

void test_invalid_dimensions()
{
    const Fixture f;
    auto d = make_direction(5);

    d.old_state.pop_back();

    bool rejected = false;

    try {
        (void)tpsf::apply_full_space_step_jvp(
            Fixture::n,
            f.old_state,
            f.new_state,
            d.old_state,
            d.new_state,
            d.old_inputs,
            d.new_inputs,
            d.reconstructed,
            f.parameters);
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    require(rejected, "Invalid tangent dimensions accepted");
}

} // namespace

int main()
{
    try {
        for (int group = 0; group <= 5; ++group)
            test_direction(group);

        test_invalid_dimensions();

        std::cout
            << "PASS full-space tangent verification\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
