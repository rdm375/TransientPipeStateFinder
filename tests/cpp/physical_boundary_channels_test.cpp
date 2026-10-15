#include "tpsf/physical_boundary_channels.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

bool close(double a, double b, double tolerance = 1e-12)
{
    return std::abs(a-b) <=
           tolerance*std::max({1.0,std::abs(a),std::abs(b)});
}

void test_extraction()
{
    constexpr int n = 3;

    const std::vector<double> state = {
        5.0e6, 12.0, 11.0, 4.9e6,
        10.0, 4.8e6, 9.0, 4.7e6
    };

    const tpsf::PhysicalBoundaryInputs inputs{
        5.0e6, 8.0
    };

    const auto channels =
        tpsf::extract_physical_boundary_channels(
            n, state, inputs);

    assert(close(channels[0], 5.0e6));
    assert(close(channels[1], 12.0));
    assert(close(channels[2], 4.7e6));
    assert(close(channels[3], 8.0));
}

void test_adjoint_identity()
{
    constexpr int n = 4;

    const std::vector<double> dx = {
        0.1, 0.2, 0.3, 0.4, 0.5,
        0.6, 0.7, 0.8, 0.9, 1.0
    };

    const tpsf::PhysicalBoundaryInputs dq{
        1.1, 1.2
    };

    const tpsf::BoundaryVector w{
        0.3, -0.7, 1.1, -1.3
    };

    const auto dh =
        tpsf::physical_boundary_channels_jvp(
            n, dx, dq);

    double lhs = 0.0;

    for (std::size_t i = 0; i < 4; ++i)
        lhs += dh[i]*w[i];

    std::vector<double> gx(2*n+2, 0.0);
    tpsf::PhysicalBoundaryInputs gq{};

    tpsf::physical_boundary_channels_vjp(
        n, w, gx, gq);

    double rhs = 0.0;

    for (std::size_t i = 0; i < dx.size(); ++i)
        rhs += dx[i]*gx[i];

    rhs += dq.inlet_pressure*gq.inlet_pressure;
    rhs += dq.outlet_flow*gq.outlet_flow;

    assert(close(lhs, rhs));
}

void test_vjp_accumulation()
{
    constexpr int n = 2;

    const tpsf::BoundaryVector w{
        1.0, 2.0, 3.0, 4.0
    };

    std::vector<double> gx(2*n+2, 7.0);

    tpsf::PhysicalBoundaryInputs gq{
        11.0, 13.0
    };

    tpsf::physical_boundary_channels_vjp(
        n, w, gx, gq);

    assert(close(gx[1], 9.0));
    assert(close(gx[2*n+1], 10.0));

    assert(close(gq.inlet_pressure, 12.0));
    assert(close(gq.outlet_flow, 17.0));

    assert(close(gx[0], 7.0));
}

void test_invalid_dimensions()
{
    std::vector<double> invalid(3, 0.0);

    bool threw = false;

    try {
        (void)tpsf::extract_physical_boundary_channels(
            3, invalid, {});
    } catch (const std::invalid_argument&) {
        threw = true;
    }

    assert(threw);
}

} // namespace

int main()
{
    test_extraction();
    test_adjoint_identity();
    test_vjp_accumulation();
    test_invalid_dimensions();
}
