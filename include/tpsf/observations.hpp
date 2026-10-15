#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace tpsf {

struct BoundaryMeasurements {
    std::vector<double> inlet_pressure;
    std::vector<double> inlet_flow;
    std::vector<double> outlet_pressure;
    std::vector<double> outlet_flow;
};

struct BoundaryMeasurementScales {
    double inlet_pressure = 1.0e4;
    double inlet_flow = 1.0;
    double outlet_pressure = 1.0e4;
    double outlet_flow = 1.0;
};

struct BoundaryObservationRecord {
    std::vector<double> time;
    BoundaryMeasurements measured;
    BoundaryMeasurementScales scales;
    bool include_initial = false;
};

inline void validate_boundary_observations(
    const BoundaryObservationRecord& record)
{
    const std::size_t n = record.time.size();

    if (n == 0)
        throw std::invalid_argument("Empty observation timeline");

    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(record.time[i]))
            throw std::invalid_argument("Nonfinite observation time");

        if (i > 0 && record.time[i] <= record.time[i - 1])
            throw std::invalid_argument(
                "Observation times must be strictly increasing");
    }

    const auto validate_channel =
        [n](const std::vector<double>& values, double scale)
    {
        if (values.size() != n)
            throw std::invalid_argument(
                "Boundary observation dimension mismatch");

        if (!std::isfinite(scale) || scale <= 0.0)
            throw std::invalid_argument(
                "Invalid measurement uncertainty");

        for (double value : values) {
            if (!std::isfinite(value))
                throw std::invalid_argument(
                    "Nonfinite boundary observation");
        }
    };

    validate_channel(
        record.measured.inlet_pressure,
        record.scales.inlet_pressure);

    validate_channel(
        record.measured.inlet_flow,
        record.scales.inlet_flow);

    validate_channel(
        record.measured.outlet_pressure,
        record.scales.outlet_pressure);

    validate_channel(
        record.measured.outlet_flow,
        record.scales.outlet_flow);
}

} // namespace tpsf
