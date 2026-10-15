#pragma once

#include "tpsf/dense_full_space_sqp.hpp"
#include "tpsf/full_space_problem.hpp"

#include <Eigen/Dense>
#include <cstddef>

namespace tpsf {

// Local, linearized observability diagnostics in scaled coordinates.
// No statistical uncertainty interpretation is implied.
struct StatefindingQuality {
    DenseSQPTermination termination =
        DenseSQPTermination::numerical_failure;

    bool converged = false;

    double objective = 0.0;
    double constraint_l2 = 0.0;
    double projected_gradient = 0.0;
    double measurement_l2 = 0.0;

    std::size_t variables = 0;
    std::size_t constraints = 0;
    std::size_t measurements = 0;

    std::size_t constraint_rank = 0;
    std::size_t feasible_dimension = 0;
    std::size_t measurement_rank = 0;
    std::size_t unresolved_directions = 0;

    // Singular values of B*Z, descending, including numerical zeros.
    Eigen::VectorXd information_singular_values;

    // Columns span the numerical nullspace of A.
    Eigen::MatrixXd feasible_basis;

    // Columns span locally unresolved feasible directions.
    Eigen::MatrixXd unresolved_basis;

    // True only if every feasible direction is locally observed.
    bool locally_identifiable = false;
};

struct StatefindingQualityOptions {
    double relative_rank_tolerance = 1.0e-12;
};

StatefindingQuality assess_statefinding_quality(
    const FullSpaceProblem& problem,
    const DenseSQPResult& solution,
    const StatefindingQualityOptions& options = {});

} // namespace tpsf
