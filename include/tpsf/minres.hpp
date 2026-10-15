#pragma once

#include <cstddef>
#include <functional>
#include <vector>

namespace tpsf {

using KrylovVector = std::vector<double>;

using SymmetricOperator =
    std::function<KrylovVector(const KrylovVector&)>;

enum class MinresStatus {
    converged,
    iteration_limit,
    breakdown,
    invalid_operator
};

struct MinresOptions {
    std::size_t max_iterations = 500;
    double relative_tolerance = 1e-10;
    double absolute_tolerance = 1e-12;
};

struct MinresResult {
    KrylovVector solution;
    MinresStatus status = MinresStatus::iteration_limit;
    std::size_t iterations = 0;
    double initial_residual_norm = 0.0;
    double final_residual_norm = 0.0;
};

MinresResult solve_minres(
    const SymmetricOperator& apply,
    const KrylovVector& rhs,
    const MinresOptions& options = {});

} // namespace tpsf
