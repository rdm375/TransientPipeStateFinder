#include "tpsf/full_space_observations.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(double a, double b, double tolerance,
          const char* message)
{
    if (!std::isfinite(a) ||
        !std::isfinite(b) ||
        std::abs(a - b) >
            tolerance * std::max({1.0, std::abs(a), std::abs(b)}))
        throw std::runtime_error(message);
}

double dot(const tpsf::BoundaryVector& a,
           const tpsf::BoundaryVector& b)
{
    double value = 0.0;

    for (std::size_t i = 0; i < tpsf::boundary_channels; ++i)
        value += a[i] * b[i];

    return value;
}

void test_full_space()
{
    tpsf::FullSpaceObservation observation{
        {5.0e6, 10.0, 4.9e6, 9.8},
        {5.0e6, 10.0, 4.9e6, 9.8},
        {1000.0, 1.0, 1000.0, 1.0},
        {5.0e6, 10.0, 4.9e6, 9.8},
        true
    };

    const auto baseline =
        tpsf::evaluate_full_space_observation(observation);

    for (double value : baseline.measurement)
        near(value, 0.0, 1.0e-14,
             "Nonzero baseline measurement residual");

    for (double value : baseline.constraint)
        near(value, 0.0, 1.0e-14,
             "Nonzero baseline constraint residual");

    // Perturb only the reconstructed inlet pressure.
    auto modified = observation;
    modified.reconstructed[0] += 100.0;

    const auto changed =
        tpsf::evaluate_full_space_observation(modified);

    near(changed.measurement[0], 0.1, 1.0e-12,
         "Incorrect measurement variation");

    near(changed.constraint[0], 100.0, 1.0e-12,
         "Incorrect constraint variation");

    // Perturb only the state inlet pressure.
    modified = observation;
    modified.state_boundary[0] += 100.0;

    const auto state_changed =
        tpsf::evaluate_full_space_observation(modified);

    near(state_changed.measurement[0], 0.0, 1.0e-12,
         "State perturbation incorrectly changed measurement");

    near(state_changed.constraint[0], -100.0, 1.0e-12,
         "Incorrect state constraint variation");

    // Jacobian/adjoint consistency.
    const tpsf::BoundaryVector delta_reconstructed{
        0.2, -0.4, 0.6, -0.8
    };

    const tpsf::BoundaryVector delta_state{
        -0.3, 0.5, -0.7, 0.9
    };

    const tpsf::FullSpaceObservationResidual dual{
        {0.3, -0.2, 0.7, -0.4},
        {-0.6, 0.8, -0.1, 0.5}
    };

    const auto tangent =
        tpsf::full_space_observation_jvp(
            observation,
            delta_reconstructed,
            delta_state);

    const auto adjoint =
        tpsf::full_space_observation_vjp(
            observation,
            dual);

    const double lhs =
        dot(tangent.measurement, dual.measurement) +
        dot(tangent.constraint, dual.constraint);

    const double rhs =
        dot(delta_reconstructed, adjoint.reconstructed) +
        dot(delta_state, adjoint.state_boundary);

    near(lhs, rhs, 1.0e-13,
         "Jacobian/adjoint dot-product identity failed");

    // Excluding an observation must not remove its
    // physical boundary constraint.
    modified = observation;
    modified.included = false;
    modified.reconstructed[0] += 100.0;

    const auto excluded =
        tpsf::evaluate_full_space_observation(modified);

    near(excluded.measurement[0], 0.0, 1.0e-14,
         "Excluded observation contributed to objective");

    near(excluded.constraint[0], 100.0, 1.0e-12,
         "Excluding observation removed physical constraint");

    std::cout
        << "PASS independent full-space observation variables\n";
}

} // namespace

int main()
{
    try {
        test_full_space();
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FAIL E2.4b.1: " << e.what() << '\n';
        return 1;
    }
}
