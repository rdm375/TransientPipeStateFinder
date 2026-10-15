#include "tpsf/full_space_adjoint.hpp"

#include "pipe_sim/physics.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace tpsf {

FullSpaceStepAdjoint apply_full_space_step_vjp(
    int n,
    std::span<const double> state_old,
    std::span<const double> state_new,
    std::span<const double> dynamics_adjoint,
    const BoundaryVector& observation_adjoint,
    const pipe_sim::TransientParameters& par)
{
    const int nu = 2*n + 2;

    if (n < 2 ||
        state_old.size() != static_cast<std::size_t>(nu) ||
        state_new.size() != static_cast<std::size_t>(nu) ||
        dynamics_adjoint.size() != static_cast<std::size_t>(nu)) {
        throw std::invalid_argument(
            "Invalid full-space adjoint dimensions");
    }

    FullSpaceStepAdjoint result;

    result.state_old.assign(nu, 0.0);
    result.state_new.assign(nu, 0.0);

    const double theta = par.theta;
    const double cold = 1.0-theta;

    const double dx = par.dx;
    const double dt = par.dt;
    const double area = par.area;
    const double diameter = par.diameter;

    const double crho = pipe_sim::density_derivative_cz(
        par.temperature, par.z, par.gas_constant);

    const double alpha = area*dx/dt;
    const double beta = 1.0/(area*dt);

    // --------------------------------------------------
    // New-state transpose:
    // F_xnew^T * lambda_F.
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
            result.state_new[col] +=
                ab[col*ldab + kv + row-col] *
                dynamics_adjoint[row];
        }
    }

    // --------------------------------------------------
    // Old-state transpose:
    // F_xold^T * lambda_F.
    // --------------------------------------------------

    // Inlet half-cell continuity.
    {
        const double lambda = dynamics_adjoint[1];

        result.state_old[0] +=
            -0.5*alpha*crho*lambda;

        result.state_old[2] += cold*lambda;
        result.state_old[1] -= cold*lambda;
    }

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

        const double lambda = dynamics_adjoint[row];

        // Momentum equation.
        result.state_old[il] +=
            cold*(-1.0/dx+dsdp)*lambda;

        result.state_old[ir] +=
            cold*(1.0/dx+dsdp)*lambda;

        result.state_old[im] +=
            (-beta+cold*dsdm)*lambda;

        // Continuity equation.
        const double lambda_c =
            dynamics_adjoint[row+1];

        if (i < n-1) {
            result.state_old[ir] +=
                -alpha*crho*lambda_c;

            result.state_old[im+2] +=
                cold*lambda_c;

            result.state_old[im] -=
                cold*lambda_c;
        } else {
            result.state_old[ir] +=
                -0.5*alpha*crho*lambda_c;

            result.state_old[im] -=
                cold*lambda_c;
        }
    }

    // --------------------------------------------------
    // Prescribed boundary input transpose.
    // --------------------------------------------------

    result.inputs_new.inlet_pressure -=
        dynamics_adjoint[0];

    result.inputs_new.outlet_flow +=
        theta*dynamics_adjoint[nu-1];

    result.inputs_old.outlet_flow +=
        cold*dynamics_adjoint[nu-1];

    // The previous prescribed inlet pressure has no
    // direct contribution to this timestep residual.

    // --------------------------------------------------
    // Observation-link transpose.
    //
    // C = b_new - H(x_new, q_new)
    //
    // Hence:
    //   g_bnew = lambda_C
    //   g_xnew -= H_xnew^T lambda_C
    //   g_qnew -= H_qnew^T lambda_C
    // --------------------------------------------------

    result.reconstructed_new = observation_adjoint;

    BoundaryVector negative_observation{};

    for (std::size_t i = 0;
         i < negative_observation.size(); ++i) {
        negative_observation[i] = -observation_adjoint[i];
    }

    physical_boundary_channels_vjp(
        n,
        negative_observation,
        result.state_new,
        result.inputs_new);

    return result;
}

} // namespace tpsf
