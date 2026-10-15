#include "tpsf/state_qualification.hpp"

#include <cassert>
#include <stdexcept>
#include <vector>

int main()
{
    using namespace tpsf;

    StateQualificationOptions options;
    options.pressure_tolerance_pa = 100.0;
    options.flow_tolerance_kg_s = 0.1;

    std::vector<StateQualificationSample> samples(5);

    const double pressure_errors[] = {
        200.0, 50.0, 150.0, 30.0, 20.0
    };

    for (std::size_t k = 0; k < samples.size(); ++k) {
        samples[k].analysis_time_s = double(k);
        samples[k].reconstruction_converged = true;
        samples[k].restart_succeeded = true;
        samples[k].pressure_error_pa = pressure_errors[k];
        samples[k].flow_error_kg_s = 0.01;
    }

    auto result =
        assess_state_qualification(samples, options);

    // A temporary acceptable state at t=1 must not
    // qualify because t=2 fails.
    assert(result.qualified);
    assert(result.qualification_time_s == 3.0);
    assert(result.samples[1].acceptable);
    assert(!result.samples[2].acceptable);

    // If the last state fails, no persistent suffix exists.
    samples.back().restart_succeeded = false;

    result = assess_state_qualification(samples, options);

    assert(!result.qualified);
    assert(!result.qualification_time_s);

    // Reject non-increasing analysis times.
    samples[2].analysis_time_s = 1.0;

    bool rejected = false;

    try {
        (void)assess_state_qualification(samples, options);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }

    assert(rejected);

    return 0;
}
