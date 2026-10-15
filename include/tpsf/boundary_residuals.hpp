#pragma once

#include "tpsf/observations.hpp"

#include <vector>

namespace tpsf {

struct BoundaryResiduals {
    std::vector<double> inlet_pressure;
    std::vector<double> inlet_flow;
    std::vector<double> outlet_pressure;
    std::vector<double> outlet_flow;
};

inline BoundaryResiduals compute_boundary_residuals(
    const BoundaryMeasurements& predicted,
    const BoundaryObservationRecord& observed)
{
    validate_boundary_observations(observed);

    const std::size_t n = observed.time.size();

    const auto check = [n](const std::vector<double>& v) {
        if (v.size() != n)
            throw std::invalid_argument(
                "Predicted boundary dimension mismatch");
    };

    check(predicted.inlet_pressure);
    check(predicted.inlet_flow);
    check(predicted.outlet_pressure);
    check(predicted.outlet_flow);

    BoundaryResiduals residuals;

    const auto compute =
        [n](const std::vector<double>& model,
            const std::vector<double>& data,
            double sigma)
    {
        std::vector<double> r(n);

        for (std::size_t i = 0; i < n; ++i) {
            if (!std::isfinite(model[i]))
                throw std::invalid_argument(
                    "Nonfinite predicted boundary value");

            r[i] = (model[i] - data[i]) / sigma;
        }

        return r;
    };

    residuals.inlet_pressure = compute(
        predicted.inlet_pressure,
        observed.measured.inlet_pressure,
        observed.scales.inlet_pressure);

    residuals.inlet_flow = compute(
        predicted.inlet_flow,
        observed.measured.inlet_flow,
        observed.scales.inlet_flow);

    residuals.outlet_pressure = compute(
        predicted.outlet_pressure,
        observed.measured.outlet_pressure,
        observed.scales.outlet_pressure);

    residuals.outlet_flow = compute(
        predicted.outlet_flow,
        observed.measured.outlet_flow,
        observed.scales.outlet_flow);

    return residuals;
}

inline double boundary_least_squares_cost(
    const BoundaryResiduals& residuals)
{
    double sum = 0.0;

    const auto accumulate = [&sum](const std::vector<double>& values) {
        for (double value : values)
            sum += value * value;
    };

    accumulate(residuals.inlet_pressure);
    accumulate(residuals.inlet_flow);
    accumulate(residuals.outlet_pressure);
    accumulate(residuals.outlet_flow);

    return 0.5 * sum;
}

} // namespace tpsf
