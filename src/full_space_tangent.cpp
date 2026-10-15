#include "tpsf/full_space_tangent.hpp"

#include "pipe_sim/physics.hpp"

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <vector>

namespace tpsf {

FullSpaceStepConstraints apply_full_space_step_jvp(
    int n,
    std::span<const double> state_old,
    std::span<const double> state_new,
    std::span<const double> tangent_old,
    std::span<const double> tangent_new,
    const PhysicalBoundaryInputs& tangent_inputs_old,
    const PhysicalBoundaryInputs& tangent_inputs_new,
    const BoundaryVector& tangent_reconstructed_new,
    const pipe_sim::TransientParameters& par)
{
    const int nu = 2*n + 2;

    if (n < 2 ||
        state_old.size() != static_cast<std::size_t>(nu) ||
        state_new.size() != static_cast<std::size_t>(nu) ||
        tangent_old.size() != static_cast<std::size_t>(nu) ||
        tangent_new.size() != static_cast<std::size_t>(nu)) {
        throw std::invalid_argument(
            "Invalid full-space tangent dimensions");
    }

    FullSpaceStepConstraints result;
    result.dynamics.assign(nu, 0.0);

    const double theta = par.theta;
    const double cold = 1.0 - theta;

    const double dx = par.dx;
    const double dt = par.dt;
    const double area = par.area;
    const double diameter = par.diameter;

    const double crho = pipe_sim::density_derivative_cz(
        par.temperature, par.z, par.gas_constant);

    const double alpha = area*dx/dt;
    const double beta = 1.0/(area*dt);

    // --------------------------------------------------
    // A * tangent_new:
    // reuse the validated analytic banded Jacobian.
    // --------------------------------------------------

    constexpr int ldab = 6;
    constexpr int kv = 3;

    std::vector<double> ab(
        static_cast<std::size_t>(ldab)*nu);

    pipe_sim::assemble_jacobian_banded(
        n, state_new, par, ab, ldab);

    for (int col = 0; col < nu; ++col) {
        const int first = std::max(0, col-1);
        const int last = std::min(nu-1, col+2);

        for (int row = first; row <= last; ++row) {
            result.dynamics[row] +=
                ab[col*ldab + kv + row-col] *
                tangent_new[col];
        }
    }

    // --------------------------------------------------
    // B * tangent_old:
    // exact previous-state derivative.
    // --------------------------------------------------

    // Inlet half-cell continuity.
    result.dynamics[1] +=
        -0.5*alpha*crho*tangent_old[0]
        + cold*(tangent_old[2]-tangent_old[1]);

    for (int i = 0; i < n; ++i) {
        const int im = 2*i+2;
        const int il = (i == 0) ? 0 : 2*i+1;
        const int ir = 2*i+3;
        const int row = 2*i+2;

        const double pleft = state_old[il];
        const double pright = state_old[ir];
        const double mdot = state_old[im];

        const double rho_left = pipe_sim::density_cz(
            pleft, par.temperature, par.z,
            par.gas_constant);

        const double rho_right = pipe_sim::density_cz(
            pright, par.temperature, par.z,
            par.gas_constant);

        const double rho_face =
            0.5*(rho_left+rho_right);

        const double re = pipe_sim::reynolds_mass(
            mdot, diameter, area, par.viscosity);

        const auto friction =
            pipe_sim::friction_sj_value_derivative(
                re, par.roughness, diameter);

        const double dfdm =
            friction.derivative_re *
            pipe_sim::reynolds_derivative_mdot(
                mdot, diameter, area, par.viscosity);

        const double sf = pipe_sim::friction_source(
            friction.value, diameter, area,
            rho_face, mdot);

        const double dsdm =
            pipe_sim::friction_source_derivative_mdot(
                friction.value, dfdm,
                diameter, area, rho_face, mdot);

        const double dsdp =
            -sf*crho/(2.0*rho_face);

        const double dp_left = tangent_old[il];
        const double dp_right = tangent_old[ir];
        const double dm = tangent_old[im];

        result.dynamics[row] +=
            -beta*dm
            + cold*(
                (dp_right-dp_left)/dx
                + dsdm*dm
                + dsdp*(dp_left+dp_right));

        if (i < n-1) {
            result.dynamics[row+1] +=
                -alpha*crho*tangent_old[ir]
                + cold*(
                    tangent_old[im+2] -
                    tangent_old[im]);
        } else {
            result.dynamics[row+1] +=
                -0.5*alpha*crho*tangent_old[ir]
                - cold*tangent_old[im];
        }
    }

    // --------------------------------------------------
    // C and D: independently prescribed boundary inputs.
    // --------------------------------------------------

    result.dynamics[0] -=
        tangent_inputs_new.inlet_pressure;

    result.dynamics[nu-1] +=
        theta*tangent_inputs_new.outlet_flow
        + cold*tangent_inputs_old.outlet_flow;

    // Previous prescribed inlet pressure does not enter
    // the timestep dynamics directly.

    // --------------------------------------------------
    // Observation links:
    // delta C = delta b - DH[delta x_new, delta q_new].
    // --------------------------------------------------

    const auto channel_tangent =
        physical_boundary_channels_jvp(
            n, tangent_new, tangent_inputs_new);

    for (std::size_t j = 0;
         j < result.observation_links.size(); ++j) {
        result.observation_links[j] =
            tangent_reconstructed_new[j]
            - channel_tangent[j];
    }

    return result;
}

} // namespace tpsf
