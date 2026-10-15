#include "tpsf/full_space_adjoint.hpp"
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

double dot(
    std::span<const double> a,
    std::span<const double> b)
{
    require(a.size() == b.size(), "Dot-product size mismatch");

    double result = 0.0;

    for (std::size_t i = 0; i < a.size(); ++i)
        result += a[i]*b[i];

    return result;
}

double boundary_dot(
    const tpsf::BoundaryVector& a,
    const tpsf::BoundaryVector& b)
{
    return dot(a, b);
}

double input_dot(
    const tpsf::PhysicalBoundaryInputs& a,
    const tpsf::PhysicalBoundaryInputs& b)
{
    return
        a.inlet_pressure*b.inlet_pressure +
        a.outlet_flow*b.outlet_flow;
}

struct Fixture {
    int n = 5;

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
};

struct Direction {
    std::vector<double> old_state;
    std::vector<double> new_state;

    tpsf::PhysicalBoundaryInputs old_inputs{};
    tpsf::PhysicalBoundaryInputs new_inputs{};

    tpsf::BoundaryVector reconstructed{};
};

Direction make_direction(int n, int group, int seed)
{
    Direction d;

    const int nu = 2*n+2;

    d.old_state.assign(nu, 0.0);
    d.new_state.assign(nu, 0.0);

    for (int i = 0; i < nu; ++i) {
        const bool pressure =
            i == 0 || (i >= 3 && i % 2 == 1);

        const double scale = pressure ? 1.0e5 : 2.0;

        if (group == 0 || group == 5) {
            d.old_state[i] =
                scale*std::sin(0.37*(i+1)*(seed+1));
        }

        if (group == 1 || group == 5) {
            d.new_state[i] =
                scale*std::cos(0.29*(i+2)*(seed+1));
        }
    }

    if (group == 2 || group == 5) {
        d.old_inputs = {
            1.0e5*std::sin(0.7*(seed+1)),
            2.0*std::cos(0.5*(seed+1))
        };
    }

    if (group == 3 || group == 5) {
        d.new_inputs = {
            1.0e5*std::cos(0.6*(seed+1)),
            2.0*std::sin(0.4*(seed+1))
        };
    }

    if (group == 4 || group == 5) {
        d.reconstructed = {
            1.0e5*std::sin(0.3*(seed+1)),
            2.0*std::cos(0.8*(seed+1)),
            1.0e5*std::cos(0.4*(seed+1)),
            2.0*std::sin(0.9*(seed+1))
        };
    }

    return d;
}

void test_transpose_identity(
    const Fixture& f,
    int group,
    int seed)
{
    const int nu = 2*f.n+2;

    const auto d = make_direction(f.n, group, seed);

    std::vector<double> lambda_f(nu);

    for (int i = 0; i < nu; ++i) {
        lambda_f[i] =
            std::sin(0.19*(i+1)*(seed+1));
    }

    tpsf::BoundaryVector lambda_c{
        std::cos(0.23*(seed+1)),
        std::sin(0.31*(seed+1)),
        std::cos(0.41*(seed+1)),
        std::sin(0.53*(seed+1))
    };

    const auto tangent =
        tpsf::apply_full_space_step_jvp(
            f.n,
            f.old_state,
            f.new_state,
            d.old_state,
            d.new_state,
            d.old_inputs,
            d.new_inputs,
            d.reconstructed,
            f.parameters);

    const auto adjoint =
        tpsf::apply_full_space_step_vjp(
            f.n,
            f.old_state,
            f.new_state,
            lambda_f,
            lambda_c,
            f.parameters);

    const double lhs =
        dot(tangent.dynamics, lambda_f) +
        boundary_dot(tangent.observation_links, lambda_c);

    const double rhs =
        dot(d.old_state, adjoint.state_old) +
        dot(d.new_state, adjoint.state_new) +
        input_dot(d.old_inputs, adjoint.inputs_old) +
        input_dot(d.new_inputs, adjoint.inputs_new) +
        boundary_dot(
            d.reconstructed,
            adjoint.reconstructed_new);

    const double scale =
        std::max({1.0, std::abs(lhs), std::abs(rhs)});

    const double relative_error =
        std::abs(lhs-rhs)/scale;

    require(std::isfinite(relative_error),
            "Nonfinite transpose error");

    require(relative_error < 1.0e-11,
            "Discrete adjoint transpose identity failed");

    std::cout
        << "group=" << group
        << " seed=" << seed
        << " relative_error=" << relative_error
        << '\n';
}

void test_invalid_dimensions()
{
    const Fixture f;

    std::vector<double> invalid_dual(
        2*f.n+1, 1.0);

    tpsf::BoundaryVector lambda_c{};

    bool rejected = false;

    try {
        (void)tpsf::apply_full_space_step_vjp(
            f.n,
            f.old_state,
            f.new_state,
            invalid_dual,
            lambda_c,
            f.parameters);
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    require(rejected,
            "Invalid adjoint dimensions accepted");
}

} // namespace

int main()
{
    try {
        const Fixture f;

        for (int group = 0; group <= 5; ++group) {
            for (int seed = 0; seed < 8; ++seed) {
                test_transpose_identity(f, group, seed);
            }
        }

        test_invalid_dimensions();

        std::cout
            << "PASS full-space discrete adjoint verification\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
