#pragma once

#include "tpsf/full_space_problem.hpp"

#include <vector>

namespace tpsf {

// Unregularized KKT operator for the equality-constrained
// linearized weighted least-squares subproblem.
//
// Variable ordering: [primal increment, constraint dual].
// All quantities use the scaled optimization coordinates.
class MatrixFreeKKT {
public:
    MatrixFreeKKT(
        const FullSpaceProblem& problem,
        std::vector<double> coordinates);

    // Apply [B^T B, A^T; A, 0] without assembling A or B.
    std::vector<double> apply(
        const std::vector<double>& vector) const;

    // Assemble -[B^T r; c], with r and c evaluated
    // at the fixed linearization point.
    std::vector<double> right_hand_side() const;

    std::size_t primal_size() const noexcept;
    std::size_t dual_size() const noexcept;

private:
    const FullSpaceProblem& problem_;
    std::vector<double> coordinates_;
};

} // namespace tpsf
