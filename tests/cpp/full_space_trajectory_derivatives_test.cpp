#include "tpsf/full_space_trajectory_derivatives.hpp"
#include "tpsf/full_space_tangent.hpp"
#include "tpsf/full_space_adjoint.hpp"
#include "tpsf/physical_boundary_channels.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace tpsf;

void require(bool ok, const std::string& message)
{
    if (!ok)
        throw std::runtime_error(message);
}

void near(double a, double b, double atol,
          double rtol, const std::string& message)
{
    require(std::isfinite(a) && std::isfinite(b),
            message + ": nonfinite");

    require(std::abs(a-b) <=
                atol + rtol *
                std::max(std::abs(a), std::abs(b)),
            message + ": disagreement: " +
            std::to_string(a) + " vs " +
            std::to_string(b));
}

template<class F>
void expect_invalid(F&& f)
{
    bool rejected = false;

    try {
        f();
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    require(rejected, "Invalid dimensions accepted");
}

void add(PhysicalBoundaryInputs& a,
         const PhysicalBoundaryInputs& b)
{
    a.inlet_pressure += b.inlet_pressure;
    a.outlet_flow += b.outlet_flow;
}

double dot(const std::vector<double>& a,
           const std::vector<double>& b)
{
    require(a.size() == b.size(), "Dot dimension mismatch");

    double result = 0.0;

    for (std::size_t i = 0; i < a.size(); ++i)
        result += a[i]*b[i];

    return result;
}

double dot(const PhysicalBoundaryInputs& a,
           const PhysicalBoundaryInputs& b)
{
    return a.inlet_pressure*b.inlet_pressure +
           a.outlet_flow*b.outlet_flow;
}

double dot(const BoundaryVector& a,
           const BoundaryVector& b)
{
    double result = 0.0;

    for (std::size_t i = 0; i < 4; ++i)
        result += a[i]*b[i];

    return result;
}

double dot(const FullSpaceTrajectoryConstraints& a,
           const FullSpaceTrajectoryConstraints& b)
{
    require(a.dynamics.size() == b.dynamics.size(),
            "Dynamics block mismatch");

    require(a.observation_links.size() ==
                b.observation_links.size(),
            "Observation block mismatch");

    double result = 0.0;

    for (std::size_t k = 0; k < a.dynamics.size(); ++k)
        result += dot(a.dynamics[k], b.dynamics[k]);

    for (std::size_t k = 0;
         k < a.observation_links.size(); ++k)
        result += dot(a.observation_links[k],
                      b.observation_links[k]);

    return result;
}

double dot(const FullSpaceTrajectoryDirection& a,
           const FullSpaceTrajectoryAdjoint& b)
{
    require(a.states.size() == b.states.size(),
            "State block mismatch");

    double result = 0.0;

    for (std::size_t k = 0; k < a.states.size(); ++k) {
        result += dot(a.states[k], b.states[k]);
        result += dot(a.inputs[k], b.inputs[k]);
        result += dot(a.reconstructed[k],
                      b.reconstructed[k]);
    }

    return result;
}

struct Fixture {
    pipe_sim::TransientParameters parameters{
        100.0, 1.0, 0.65,
        0.5, 0.19634954084936207,
        288.15, 0.9, 518.3,
        1.1e-5, 1.0e-5
    };

    FullSpaceTrajectory trajectory;
    FullSpaceTrajectoryDirection direction;
    FullSpaceTrajectoryConstraints dual;

    Fixture(int n, int steps)
    {
        trajectory.n = n;

        const int nu = 2*n+2;

        for (int k = 0; k <= steps; ++k) {
            std::vector<double> x(nu);
            std::vector<double> dx(nu);

            for (int i = 0; i < nu; ++i) {
                if (i == 0 || (i >= 3 && i%2 == 1)) {
                    x[i] = 5.0e6 - 2000.0*i +
                           500.0*k;

                    dx[i] = 100.0 *
                            std::sin(0.37*(i+1)*(k+1));
                } else {
                    x[i] = 10.0 +
                           0.08*i +
                           0.03*k;

                    dx[i] = 0.1 *
                            std::cos(0.29*(i+1)*(k+1));
                }
            }

            trajectory.states.push_back(x);
            direction.states.push_back(dx);

            const PhysicalBoundaryInputs q{
                5.0e6 + 700.0*k,
                10.0 + 0.13*k
            };

            const PhysicalBoundaryInputs dq{
                100.0*std::cos(0.4*(k+1)),
                0.1*std::sin(0.6*(k+1))
            };

            trajectory.inputs.push_back(q);
            direction.inputs.push_back(dq);

            const auto physical =
                extract_physical_boundary_channels(
                    n, x, q);

            BoundaryVector b = physical;
            BoundaryVector db{};

            for (std::size_t c = 0; c < 4; ++c) {
                b[c] += 0.05*(k+1)*(c+1);

                // Pressure channels have pressure-scale
                // perturbations; flow channels use flow scale.
                db[c] = (c == 0 || c == 2 ? 100.0 : 0.1)
                      * std::sin(0.2*(k+1)*(c+1));
            }

            trajectory.reconstructed.push_back(b);
            direction.reconstructed.push_back(db);
        }

        for (int k = 0; k < steps; ++k) {
            std::vector<double> lambda(nu);

            for (int i = 0; i < nu; ++i)
                lambda[i] =
                    std::sin(0.17*(k+1)*(i+1));

            dual.dynamics.push_back(lambda);
        }

        for (int k = 0; k <= steps; ++k) {
            BoundaryVector mu{};

            for (std::size_t c = 0; c < 4; ++c)
                mu[c] =
                    std::cos(0.31*(k+1)*(c+1));

            dual.observation_links.push_back(mu);
        }
    }
};

FullSpaceTrajectory perturb(
    const FullSpaceTrajectory& original,
    const FullSpaceTrajectoryDirection& direction,
    double scale)
{
    auto result = original;

    for (std::size_t k = 0; k < result.states.size(); ++k) {
        for (std::size_t i = 0;
             i < result.states[k].size(); ++i) {
            result.states[k][i] +=
                scale*direction.states[k][i];
        }

        result.inputs[k].inlet_pressure +=
            scale*direction.inputs[k].inlet_pressure;

        result.inputs[k].outlet_flow +=
            scale*direction.inputs[k].outlet_flow;

        for (std::size_t c = 0; c < 4; ++c) {
            result.reconstructed[k][c] +=
                scale*direction.reconstructed[k][c];
        }
    }

    return result;
}

void test_transpose(int n, int steps)
{
    const Fixture f(n, steps);

    const auto tangent =
        apply_full_space_trajectory_jvp(
            f.trajectory, f.direction, f.parameters);

    const auto adjoint =
        apply_full_space_trajectory_vjp(
            f.trajectory, f.dual, f.parameters);

    const double lhs = dot(tangent, f.dual);
    const double rhs = dot(f.direction, adjoint);

    near(lhs, rhs, 1.0e-9, 2.0e-12,
         "Trajectory transpose identity");
}

void test_step_assembly()
{
    const Fixture f(5, 4);

    const auto tangent =
        apply_full_space_trajectory_jvp(
            f.trajectory, f.direction, f.parameters);

    for (std::size_t k = 0;
         k < f.trajectory.steps(); ++k) {

        const auto step =
            apply_full_space_step_jvp(
                f.trajectory.n,
                f.trajectory.states[k],
                f.trajectory.states[k+1],
                f.direction.states[k],
                f.direction.states[k+1],
                f.direction.inputs[k],
                f.direction.inputs[k+1],
                f.direction.reconstructed[k+1],
                f.parameters);

        require(tangent.dynamics[k] == step.dynamics,
                "Trajectory dynamics differ from step JVP");

        require(tangent.observation_links[k+1] ==
                    step.observation_links,
                "Trajectory observation differs from step JVP");
    }

    const auto dh0 =
        physical_boundary_channels_jvp(
            f.trajectory.n,
            f.direction.states[0],
            f.direction.inputs[0]);

    for (std::size_t c = 0; c < 4; ++c) {
        near(tangent.observation_links[0][c],
             f.direction.reconstructed[0][c]-dh0[c],
             0.0, 0.0,
             "Initial observation tangent");
    }
}

void test_adjoint_accumulation()
{
    const Fixture f(5, 4);

    const auto adjoint =
        apply_full_space_trajectory_vjp(
            f.trajectory, f.dual, f.parameters);

    FullSpaceTrajectoryAdjoint expected;

    const auto count = f.trajectory.states.size();
    const auto nu =
        static_cast<std::size_t>(2*f.trajectory.n+2);

    expected.states.assign(count,
                           std::vector<double>(nu, 0.0));
    expected.inputs.resize(count);
    expected.reconstructed.resize(count);

    expected.reconstructed[0] =
        f.dual.observation_links[0];

    BoundaryVector negative_initial{};

    for (std::size_t c = 0; c < 4; ++c)
        negative_initial[c] =
            -f.dual.observation_links[0][c];

    physical_boundary_channels_vjp(
        f.trajectory.n,
        negative_initial,
        expected.states[0],
        expected.inputs[0]);

    for (std::size_t k = 0;
         k < f.trajectory.steps(); ++k) {

        const auto step =
            apply_full_space_step_vjp(
                f.trajectory.n,
                f.trajectory.states[k],
                f.trajectory.states[k+1],
                f.dual.dynamics[k],
                f.dual.observation_links[k+1],
                f.parameters);

        for (std::size_t i = 0; i < nu; ++i) {
            expected.states[k][i] += step.state_old[i];
            expected.states[k+1][i] += step.state_new[i];
        }

        add(expected.inputs[k], step.inputs_old);
        add(expected.inputs[k+1], step.inputs_new);

        expected.reconstructed[k+1] =
            step.reconstructed_new;
    }

    for (std::size_t k = 0; k < count; ++k) {
        require(adjoint.states[k] == expected.states[k],
                "Shared-state adjoint accumulation incorrect");

        near(adjoint.inputs[k].inlet_pressure,
             expected.inputs[k].inlet_pressure,
             0.0, 0.0,
             "Inlet-pressure adjoint accumulation");

        near(adjoint.inputs[k].outlet_flow,
             expected.inputs[k].outlet_flow,
             0.0, 0.0,
             "Outlet-flow adjoint accumulation");

        require(adjoint.reconstructed[k] ==
                    expected.reconstructed[k],
                "Reconstructed adjoint accumulation");
    }
}

void test_finite_difference()
{
    const Fixture f(5, 3);

    const auto tangent =
        apply_full_space_trajectory_jvp(
            f.trajectory, f.direction, f.parameters);

    const auto analytic = tangent.combined();

    // Physically scaled direction:
    // pressure ~100 Pa, mass flow ~0.1 kg/s.
    //
    // h=0.01 gives pressure perturbations ~1 Pa,
    // avoiding excessively small differences relative
    // to the 5 MPa base pressure.
    constexpr double h = 0.01;

    const auto plus = perturb(
        f.trajectory, f.direction, h);

    const auto minus = perturb(
        f.trajectory, f.direction, -h);

    const auto gp =
        assemble_full_space_trajectory_constraints(
            plus, f.parameters).combined();

    const auto gm =
        assemble_full_space_trajectory_constraints(
            minus, f.parameters).combined();

    require(gp.size() == analytic.size(),
            "Finite-difference dimension mismatch");

    double maximum_relative_error = 0.0;

    for (std::size_t i = 0; i < gp.size(); ++i) {
        const double numerical =
            (gp[i]-gm[i])/(2.0*h);

        const double scale =
            std::max({1.0,
                      std::abs(numerical),
                      std::abs(analytic[i])});

        const double relative_error =
            std::abs(numerical-analytic[i])/scale;

        maximum_relative_error =
            std::max(maximum_relative_error,
                     relative_error);
    }

    std::cout
        << "Finite-difference maximum scaled error: "
        << maximum_relative_error << '\n';

    require(maximum_relative_error < 2.0e-5,
            "Trajectory finite-difference verification failed");
}

void test_invalid_dimensions()
{
    const Fixture f(5, 3);

    auto bad_direction = f.direction;
    bad_direction.states[1].pop_back();

    expect_invalid([&] {
        apply_full_space_trajectory_jvp(
            f.trajectory, bad_direction, f.parameters);
    });

    auto bad_dual = f.dual;
    bad_dual.dynamics[0].pop_back();

    expect_invalid([&] {
        apply_full_space_trajectory_vjp(
            f.trajectory, bad_dual, f.parameters);
    });

    auto bad_trajectory = f.trajectory;
    bad_trajectory.inputs.pop_back();

    expect_invalid([&] {
        apply_full_space_trajectory_jvp(
            bad_trajectory, f.direction, f.parameters);
    });
}

} // namespace

int main()
{
    try {
        for (int n : {2, 5, 10}) {
            for (int steps : {1, 2, 4, 8}) {
                test_transpose(n, steps);
            }
        }

        test_step_assembly();
        test_adjoint_accumulation();
        test_finite_difference();
        test_invalid_dimensions();

        std::cout
            << "PASS E2.4b.8 trajectory derivatives\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E2.4b.8: "
            << error.what() << '\n';

        return 1;
    }
}
