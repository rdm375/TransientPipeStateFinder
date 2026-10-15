#include "tpsf/scaled_initial_state.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void check(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void check_close(double actual, double expected, double tolerance)
{
    if (!std::isfinite(actual) ||
        std::abs(actual - expected) > tolerance)
        throw std::runtime_error("Numerical comparison failed");
}

} // namespace

int main()
{
    try {
        using namespace tpsf;

        const InitialState baseline{
            {8.0e6, 7.9e6, 7.8e6},
            {100.0, 100.0},
            100.0
        };

        const InitialState state{
            {8.0e6, 7.95e6, 7.7e6},
            {105.0, 95.0},
            110.0
        };

        // Forward coordinate transformation.
        const auto x = pack_scaled_state(state, baseline);

        const std::vector<double> expected{
            0.5, -1.0, 0.5, -0.5, 1.0
        };

        check(x.size() == expected.size(),
              "Incorrect scaled-state dimension");

        for (std::size_t i = 0; i < x.size(); ++i)
            check_close(x[i], expected[i], 1.0e-12);

        // Inverse transformation.
        const auto recovered = unpack_scaled_state(x, baseline);

        for (std::size_t i = 0; i < state.pressure.size(); ++i)
            check_close(recovered.pressure[i],
                        state.pressure[i], 1.0e-8);

        for (std::size_t i = 0; i < state.flow.size(); ++i)
            check_close(recovered.flow[i],
                        state.flow[i], 1.0e-12);

        check_close(recovered.inlet_flow,
                    state.inlet_flow, 1.0e-12);

        // Prescribed inlet pressure must remain unchanged.
        check_close(recovered.pressure[0],
                    baseline.pressure[0], 0.0);

        // Chain-rule gradient transformation.
        const std::vector<double> gradient{
            1.0, 2.0, 3.0, 4.0, 5.0
        };

        const auto scaled =
            scale_adjoint_gradient(gradient, 2);

        const std::vector<double> expected_gradient{
            1.0e5, 2.0e5, 30.0, 40.0, 50.0
        };

        for (std::size_t i = 0; i < scaled.size(); ++i)
            check_close(scaled[i],
                        expected_gradient[i], 1.0e-10);

        // Reject inconsistent dimensions.
        bool rejected = false;

        try {
            unpack_scaled_state({1.0}, baseline);
        }
        catch (const std::invalid_argument&) {
            rejected = true;
        }

        check(rejected, "Invalid dimensions were accepted");

        // Reject invalid scaling.
        rejected = false;

        try {
            pack_scaled_state(
                state, baseline, InitialStateScales{0.0, 10.0});
        }
        catch (const std::invalid_argument&) {
            rejected = true;
        }

        check(rejected, "Invalid pressure scale was accepted");

        std::cout
            << "PASS scaled initial-state representation\n";

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
