#pragma once

#include "tpsf/physical_boundary_channels.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace tpsf {

struct StateQualityMetrics {
    // Numerical consistency.
    std::optional<double> dynamics_residual_linf;
    std::optional<double> mass_conservation_defect;

    // Physical operating diagnostics.
    std::optional<double> mass_flow_imbalance;
    std::optional<double> linepack;

    // Measurement agreement.
    std::optional<double> measurement_weighted_rms;
    std::optional<double> measurement_max_normalized_error;
    std::size_t measurement_count = 0;

    // Reconstruction uncertainty.
    std::optional<double> state_uncertainty_rms;
    std::optional<double> state_uncertainty_max;
};

struct ScoredStateSnapshot {
    std::size_t step = 0;
    double time = 0.0;

    // [p_0, m_in, m_0, p_1, ..., m_(n-1), p_n]
    std::vector<double> state;

    PhysicalBoundaryInputs inputs{};
    BoundaryVector physical_channels{};

    StateQualityMetrics quality;
};

struct ScoredTrajectory {
    int cells = 0;
    std::vector<ScoredStateSnapshot> snapshots;
};

} // namespace tpsf
