#pragma once

#include "tpsf/statefinding_quality.hpp"

#include <Eigen/Dense>
#include <cstddef>
#include <vector>

namespace tpsf {

enum class EngineeringQuantity {
    initial_linepack,
    terminal_linepack,
    minimum_pressure
};

struct QuantitySensitivity {
    EngineeringQuantity quantity;

    // kg for linepack; Pa for pressure.
    double value = 0.0;

    // Derivatives per unit displacement in scaled coordinates.
    Eigen::VectorXd gradient;
    Eigen::VectorXd feasible_sensitivity;
    Eigen::VectorXd unresolved_sensitivity;

    double feasible_sensitivity_norm = 0.0;
    double unresolved_sensitivity_norm = 0.0;

    // Minimum pressure may have multiple active locations.
    // In that case no unique gradient is reported.
    bool differentiable = true;
    std::size_t active_minimum_count = 0;
};

struct EngineeringSensitivityOptions {
    // Dimensionless displacement in optimization coordinates.
    double difference_step = 1.0e-5;

    // Absolute pressure tolerance for identifying tied minima.
    double minimum_pressure_tie_tolerance_pa = 1.0e-6;
};

double evaluate_engineering_quantity(
    const FullSpaceProblem& problem,
    const pipe_sim::TransientParameters& parameters,
    const Eigen::VectorXd& coordinates,
    EngineeringQuantity quantity);

QuantitySensitivity assess_engineering_sensitivity(
    const FullSpaceProblem& problem,
    const pipe_sim::TransientParameters& parameters,
    const Eigen::VectorXd& coordinates,
    const StatefindingQuality& quality,
    EngineeringQuantity quantity,
    const EngineeringSensitivityOptions& options = {});

} // namespace tpsf
