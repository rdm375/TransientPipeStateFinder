// E4.2b: dense diagnostic oracle for Jacobian structure and rank.
//
// Reuse the existing synthetic configuration and its established
// variable ordering. This executable is not a production solver.

#define main tpsf_existing_dense_test_main
#include "dense_full_space_subproblem_test.cpp"
#undef main

#include <Eigen/SVD>

#include <iomanip>
#include <limits>
#include <string>

namespace {

struct RankReport {
    Eigen::Index rank = 0;
    double largest = 0.0;
    double smallest_retained = 0.0;
    double threshold = 0.0;
};

RankReport numerical_rank(
    const Eigen::MatrixXd& matrix,
    double relative_tolerance)
{
    RankReport result;

    if (matrix.rows() == 0 || matrix.cols() == 0)
        return result;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);

    const auto singular = svd.singularValues();

    if (singular.size() == 0)
        return result;

    result.largest = singular[0];
    result.threshold =
        relative_tolerance * result.largest;

    for (Eigen::Index i = 0; i < singular.size(); ++i) {
        if (singular[i] > result.threshold) {
            ++result.rank;
            result.smallest_retained = singular[i];
        }
    }

    return result;
}

void print_rank(
    const std::string& name,
    const RankReport& report)
{
    std::cout
        << name
        << " rank=" << report.rank
        << " sigma_max=" << report.largest
        << " sigma_min_retained="
        << report.smallest_retained
        << " threshold=" << report.threshold;

    if (report.rank > 0)
        std::cout
            << " retained_condition_estimate="
            << report.largest/report.smallest_retained;

    std::cout << '\n';
}


void diagnose_rank_sensitivity(
    const std::string& name,
    const tpsf::FullSpaceProblem& problem)
{
    using namespace tpsf;

    const Eigen::Index n =
        static_cast<Eigen::Index>(
            problem.free_variable_count());

    const Eigen::Index m =
        static_cast<Eigen::Index>(
            problem.constraint_count());

    const auto dense =
        assemble_dense_full_space_subproblem(
            problem, Eigen::VectorXd::Zero(n));

    const Eigen::MatrixXd& A =
        dense.constraint_jacobian;

    const Eigen::MatrixXd& B =
        dense.measurement_jacobian;

    Eigen::MatrixXd C(A.rows()+B.rows(), n);

    C.topRows(A.rows()) = A;
    C.bottomRows(B.rows()) = B;

    Eigen::MatrixXd K =
        Eigen::MatrixXd::Zero(n+m, n+m);

    K.topLeftCorner(n,n) = B.transpose()*B;
    K.topRightCorner(n,m) = A.transpose();
    K.bottomLeftCorner(m,n) = A;

    std::cout
        << "\n=== E4.2c rank sensitivity: "
        << name << " ===\n";

    constexpr double tolerances[] = {
        1e-6, 1e-8, 1e-10, 1e-12, 1e-14
    };

    std::cout
        << "tolerance rank_A rank_C rank_K "
           "expected_K_nullity observed_K_nullity\n";

    for (double tolerance : tolerances) {
        const auto ra = numerical_rank(A, tolerance);
        const auto rc = numerical_rank(C, tolerance);
        const auto rk = numerical_rank(K, tolerance);

        const Eigen::Index expected_nullity =
            n-rc.rank + m-ra.rank;

        const Eigen::Index observed_nullity =
            n+m-rk.rank;

        std::cout
            << tolerance << ' '
            << ra.rank << ' '
            << rc.rank << ' '
            << rk.rank << ' '
            << expected_nullity << ' '
            << observed_nullity << '\n';
    }

    if (n == 0) {
        std::cout
            << "nullspace_analysis_skipped: "
               "no free variables\n";
        return;
    }

    // Full V is required to recover the null space of
    // a rectangular constraint Jacobian.
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        A, Eigen::ComputeFullV);

    const Eigen::VectorXd singular =
        svd.singularValues();

    const double sigma_max =
        singular.size() ? singular[0] : 0.0;

    std::cout << "A singular values:\n";

    for (Eigen::Index i = 0; i < singular.size(); ++i)
        std::cout << i << ' ' << singular[i] << '\n';

    constexpr double null_tolerance = 1e-10;

    const auto ra =
        numerical_rank(A, null_tolerance);

    const Eigen::Index nullity =
        n-ra.rank;

    const Eigen::MatrixXd Z =
        svd.matrixV().rightCols(nullity);

    const Eigen::MatrixXd AZ = A*Z;
    const Eigen::MatrixXd BZ = B*Z;

    std::cout
        << "numerical_nullity_A=" << nullity
        << " norm_AZ=" << AZ.norm()
        << " norm_BZ=" << BZ.norm()
        << '\n';

    const auto rbz =
        numerical_rank(BZ, null_tolerance);

    print_rank("BZ", rbz);

    std::cout
        << "BZ_identified_directions="
        << rbz.rank
        << " BZ_unidentified_directions="
        << nullity-rbz.rank
        << '\n';

    if (BZ.rows() > 0 && BZ.cols() > 0) {
        Eigen::JacobiSVD<Eigen::MatrixXd> bsvd(BZ);

        std::cout << "BZ singular values:\n";

        const auto values = bsvd.singularValues();

        for (Eigen::Index i = 0; i < values.size(); ++i)
            std::cout << i << ' ' << values[i] << '\n';
    }
}

void diagnose(
    const std::string& name,
    const tpsf::FullSpaceProblem& problem,
    int cells,
    int steps,
    double tolerance)
{
    using namespace tpsf;

    const auto n =
        static_cast<Eigen::Index>(
            problem.free_variable_count());

    const auto m =
        static_cast<Eigen::Index>(
            problem.constraint_count());

    const Eigen::VectorXd u =
        Eigen::VectorXd::Zero(n);

    const auto dense =
        assemble_dense_full_space_subproblem(problem, u);

    const Eigen::MatrixXd& A =
        dense.constraint_jacobian;

    const Eigen::MatrixXd& B =
        dense.measurement_jacobian;

    Eigen::MatrixXd stacked(A.rows()+B.rows(), n);

    stacked.topRows(A.rows()) = A;
    stacked.bottomRows(B.rows()) = B;

    const auto rank_A =
        numerical_rank(A, tolerance);

    const auto rank_B =
        numerical_rank(B, tolerance);

    const auto rank_stacked =
        numerical_rank(stacked, tolerance);

    const Eigen::Index feasible_dimension =
        n-rank_A.rank;

    const Eigen::Index unidentified_dimension =
        n-rank_stacked.rank;

    const Eigen::Index identified_dimension =
        rank_stacked.rank-rank_A.rank;

    // KKT rank is reported independently. It can expose
    // redundant constraint multipliers even when the primal
    // state is locally identifiable.
    Eigen::MatrixXd K =
        Eigen::MatrixXd::Zero(n+m, n+m);

    K.topLeftCorner(n, n) =
        B.transpose()*B;

    K.topRightCorner(n, m) = A.transpose();
    K.bottomLeftCorner(m, n) = A;

    const auto rank_K =
        numerical_rank(K, tolerance);

    std::cout
        << "\n=== " << name << " ===\n"
        << "cells=" << cells
        << " steps=" << steps
        << " free_variables=" << n
        << " constraints=" << m
        << " measurements=" << B.rows()
        << '\n';

    print_rank("A", rank_A);
    print_rank("B", rank_B);
    print_rank("[A;B]", rank_stacked);
    print_rank("KKT", rank_K);

    std::cout
        << "feasible_nullity=" << feasible_dimension
        << " identified_feasible_directions="
        << identified_dimension
        << " unidentified_feasible_directions="
        << unidentified_dimension
        << " redundant_constraint_rows="
        << m-rank_A.rank
        << " KKT_nullity=" << n+m-rank_K.rank
        << '\n';

    // Structural coupling study. A time-level block consists
    // of state, prescribed inputs, and reconstructed channels.
    // Free-variable columns retain their original time level.
    //
    // The fixture fixes inlet pressure at each level.
    const int width = 2*cells+2;
    const int block = width+6;
    const int levels = steps+1;

    std::vector<int> column_level;

    for (int k = 0; k < levels; ++k) {
        for (int j = 0; j < block; ++j) {
            if (j == width)
                continue;

            column_level.push_back(k);
        }
    }

    if (static_cast<Eigen::Index>(column_level.size()) != n) {
        std::cout
            << "structure_skipped: fixture free-variable "
               "layout differs from expected ordering\n";
        return;
    }

    const double structural_threshold =
        1.0e-12 *
        std::max(1.0, A.cwiseAbs().maxCoeff());

    Eigen::MatrixXi dynamics_blocks =
        Eigen::MatrixXi::Zero(steps, levels);

    Eigen::MatrixXi observation_blocks =
        Eigen::MatrixXi::Zero(levels, levels);

    // Dynamics occupy the first steps*width rows.
    // Observation links follow, four rows per time level.
    for (Eigen::Index row = 0; row < A.rows(); ++row) {
        const bool dynamics =
            row < static_cast<Eigen::Index>(steps*width);

        const int row_level = dynamics
            ? static_cast<int>(row/width)
            : static_cast<int>(
                (row-steps*width)/4);

        for (Eigen::Index col = 0; col < A.cols(); ++col) {
            if (std::abs(A(row,col)) <= structural_threshold)
                continue;

            const int col_level =
                column_level[static_cast<std::size_t>(col)];

            if (dynamics)
                ++dynamics_blocks(row_level,col_level);
            else
                ++observation_blocks(row_level,col_level);
        }
    }

    std::cout
        << "structural_threshold="
        << structural_threshold << '\n';

    std::cout << "dynamics_time_blocks:\n"
              << dynamics_blocks << '\n';

    std::cout << "observation_time_blocks:\n"
              << observation_blocks << '\n';

    // Check the expected adjacent-level dynamics coupling
    // and same-level observation-link coupling.
    int off_band = 0;

    for (int k = 0; k < steps; ++k)
        for (int j = 0; j < levels; ++j)
            if (j != k && j != k+1)
                off_band += dynamics_blocks(k,j);

    for (int k = 0; k < levels; ++k)
        for (int j = 0; j < levels; ++j)
            if (j != k)
                off_band += observation_blocks(k,j);

    std::cout
        << "off_band_entries=" << off_band << '\n';

    if (off_band != 0)
        throw std::runtime_error(
            "Unexpected off-band Jacobian coupling");
}

} // namespace

int main()
{
    try {
        constexpr double tolerance = 1.0e-10;

        std::cout << std::scientific
                  << std::setprecision(12);

        const tpsf::FullSpaceProblem ordinary(
            make_configuration());

        const tpsf::FullSpaceProblem unobserved(
            make_configuration(false, true));

        const tpsf::FullSpaceProblem fixed(
            make_configuration(true, false));

        diagnose(
            "ordinary", ordinary, 2, 2, tolerance);

        diagnose(
            "no measurements", unobserved, 2, 2,
            tolerance);

        diagnose(
            "all fixed", fixed, 2, 2, tolerance);

        diagnose_rank_sensitivity(
            "ordinary", ordinary);

        diagnose_rank_sensitivity(
            "no measurements", unobserved);

        std::cout
            << "\nPASS E4.2b/E4.2c KKT structure diagnostic\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E4.2b: "
            << error.what() << '\n';

        return 1;
    }
}
