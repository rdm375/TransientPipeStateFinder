#include "tpsf/dense_full_space_sqp.hpp"
#include "tpsf/dense_full_space_subproblem.hpp"
#include "tpsf/scheduled_forward.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace tpsf;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

double maximum_difference(
    const std::vector<double>& a,
    const std::vector<double>& b)
{
    require(a.size() == b.size(), "Comparison size mismatch");

    double maximum = 0.0;

    for (std::size_t i = 0; i < a.size(); ++i)
        maximum = std::max(maximum, std::abs(a[i]-b[i]));

    return maximum;
}

Configuration make_configuration()
{
    constexpr int n = 2;
    constexpr double diameter = 0.5;
    constexpr double area = 0.19634954084936207;

    return {
        n,
        2,
        pipe_sim::TransientParameters{
            100.0,
            1.0,
            0.65,
            diameter,
            area,
            288.15,
            0.9,
            518.3,
            1.1e-5,
            1.0e-5
        },
        pipe_sim::DemandRamp{
            5.0e6, 10.0, 10.2, 2.0
        },
        pipe_sim::NewtonOptions{
            1.0e-8, 1.0e-10, 30, false
        }
    };
}

FullSpaceTrajectory make_reference(
    const Configuration& cfg,
    const RecordedForwardResult& forward)
{
    FullSpaceTrajectory truth;
    truth.n = cfg.cells;

    for (const auto& snapshot : forward.trajectory.snapshots) {
        truth.states.push_back(snapshot.state);
        truth.inputs.push_back(snapshot.inputs);
        truth.reconstructed.push_back(snapshot.physical_channels);
    }

    require(
        truth.states.size() ==
            static_cast<std::size_t>(cfg.steps+1),
        "Incomplete recorded trajectory");

    return truth;
}

double maximum_constraint_residual(
    const FullSpaceTrajectory& trajectory,
    const pipe_sim::TransientParameters& parameters)
{
    const auto constraints =
        assemble_full_space_trajectory_constraints(
            trajectory, parameters);

    const auto combined = constraints.combined();

    double maximum = 0.0;

    for (double value : combined)
        maximum = std::max(maximum, std::abs(value));

    return maximum;
}


struct ConservationAudit {
    double maximum_mass_imbalance = 0.0;
    double maximum_linepack_change_rate = 0.0;
    double maximum_conservation_defect = 0.0;
    double maximum_integrated_defect = 0.0;
};

ConservationAudit audit_conservation(
    const FullSpaceTrajectory& trajectory,
    const pipe_sim::TransientParameters& parameters)
{
    ConservationAudit audit;

    const int n = trajectory.n;
    const double dt = parameters.dt;
    const double theta = parameters.theta;

    std::vector<double> linepack;
    std::vector<double> imbalance;

    linepack.reserve(trajectory.states.size());
    imbalance.reserve(trajectory.states.size());

    for (std::size_t k = 0;
         k < trajectory.states.size(); ++k) {

        const auto& state = trajectory.states[k];
        const auto& inputs = trajectory.inputs[k];

        std::vector<double> pressure(n+1);

        pressure[0] = state[0];

        for (int i = 0; i < n; ++i)
            pressure[i+1] = state[2*i+3];

        const double mass =
            pipe_sim::linepack_cz(
                n, pressure, parameters);

        const double mass_imbalance =
            state[1] - inputs.outlet_flow;

        linepack.push_back(mass);
        imbalance.push_back(mass_imbalance);

        audit.maximum_mass_imbalance =
            std::max(
                audit.maximum_mass_imbalance,
                std::abs(mass_imbalance));
    }

    double integrated_boundary_exchange = 0.0;

    for (std::size_t k = 1;
         k < linepack.size(); ++k) {

        const double inventory_rate =
            (linepack[k]-linepack[k-1])/dt;

        const double boundary_rate =
            theta*imbalance[k] +
            (1.0-theta)*imbalance[k-1];

        const double defect =
            inventory_rate-boundary_rate;

        integrated_boundary_exchange +=
            dt*boundary_rate;

        const double integrated_defect =
            (linepack[k]-linepack[0]) -
            integrated_boundary_exchange;

        audit.maximum_linepack_change_rate =
            std::max(
                audit.maximum_linepack_change_rate,
                std::abs(inventory_rate));

        audit.maximum_conservation_defect =
            std::max(
                audit.maximum_conservation_defect,
                std::abs(defect));

        audit.maximum_integrated_defect =
            std::max(
                audit.maximum_integrated_defect,
                std::abs(integrated_defect));
    }

    return audit;
}


struct RankAudit {
    Eigen::Index variables = 0;
    Eigen::Index constraints = 0;
    Eigen::Index measurements = 0;

    Eigen::Index constraint_rank = 0;
    Eigen::Index nullity = 0;
    Eigen::Index reduced_measurement_rank = 0;
    Eigen::Index unresolved = 0;
};

Eigen::Index numerical_rank(
    const Eigen::VectorXd& singular_values,
    double relative_tolerance)
{
    if (singular_values.size() == 0)
        return 0;

    const double largest = singular_values[0];

    if (largest == 0.0)
        return 0;

    Eigen::Index rank = 0;

    for (Eigen::Index i = 0;
         i < singular_values.size(); ++i) {

        if (singular_values[i] >
            relative_tolerance*largest)
            ++rank;
    }

    return rank;
}

RankAudit audit_rank(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& coordinates)
{
    constexpr double relative_rank_tolerance = 1.0e-12;

    const auto subproblem =
        assemble_dense_full_space_subproblem(
            problem, coordinates);

    const auto& A = subproblem.constraint_jacobian;
    const auto& B = subproblem.measurement_jacobian;

    RankAudit audit;

    audit.variables = A.cols();
    audit.constraints = A.rows();
    audit.measurements = B.rows();

    const Eigen::Index n = A.cols();

    if (n == 0)
        return audit;

    // Full right singular vectors are required to construct
    // a complete basis for the constraint null space.
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        A, Eigen::ComputeFullV);

    audit.constraint_rank =
        numerical_rank(
            svd.singularValues(),
            relative_rank_tolerance);

    audit.nullity =
        n-audit.constraint_rank;

    if (audit.nullity == 0)
        return audit;

    const Eigen::MatrixXd Z =
        svd.matrixV().rightCols(audit.nullity);

    const Eigen::MatrixXd reduced = B*Z;

    const Eigen::JacobiSVD<Eigen::MatrixXd> reduced_svd(
        reduced, Eigen::ComputeThinU |
                 Eigen::ComputeThinV);

    audit.reduced_measurement_rank =
        numerical_rank(
            reduced_svd.singularValues(),
            relative_rank_tolerance);

    audit.unresolved =
        audit.nullity-audit.reduced_measurement_rank;

    return audit;
}

void test_synthetic_reconstruction()
{
    const auto cfg = make_configuration();

    InitialState initial{
        std::vector<double>(cfg.cells+1, 5.0e6),
        std::vector<double>(cfg.cells, 10.0),
        10.0
    };

    BoundarySchedule schedule{
        {5.0e6, 5.0e6, 5.0e6},
        {10.0, 10.1, 10.2}
    };

    const auto forward =
        run_forward_scheduled_recorded(
            cfg, initial, schedule);

    require(forward.summary.info == 0,
            "Forward solver failed");

    require(forward.summary.completed_steps == cfg.steps,
            "Forward solver did not complete");

    const auto truth = make_reference(cfg, forward);

    const double reference_defect =
        maximum_constraint_residual(
            truth, cfg.parameters);

    std::cout << std::scientific
              << std::setprecision(12);

    std::cout
        << "Reference constraint Linf: "
        << reference_defect << '\n';

    // This is a diagnostic gate, not a tolerance adjustment.
    // The forward integrator must satisfy the exact constraint
    // assembler before reconstruction is attempted.
    require(reference_defect < 1.0e-5,
            "Forward trajectory and constraint assembler disagree");

    FullSpaceProblemConfiguration problem_cfg;
    problem_cfg.reference = truth;
    problem_cfg.parameters = cfg.parameters;

    const std::size_t levels = truth.states.size();
    const std::size_t state_size =
        static_cast<std::size_t>(2*cfg.cells+2);

    const std::size_t variables_per_level =
        state_size+2+4;

    const std::size_t physical_variables =
        levels*variables_per_level;

    problem_cfg.variables.assign(
        physical_variables, ProblemVariable{});

    // Fix the known initial state and all prescribed inputs.
    for (std::size_t k = 0; k < levels; ++k) {
        const std::size_t offset =
            k*variables_per_level;

        for (std::size_t j = 0; j < state_size; ++j) {
            auto& variable =
                problem_cfg.variables[offset+j];

            variable.free = (k != 0);

            // Pressure and mass flow have different scales.
            variable.scale =
                (j == 0 || (j >= 3 && j % 2 == 1))
                    ? 1.0e6
                    : 10.0;
        }

        // Prescribed inlet pressure and outlet mass flow.
        problem_cfg.variables[offset+state_size].free = false;
        problem_cfg.variables[offset+state_size+1].free = false;

        // Reconstructed physical boundary channels.
        for (std::size_t channel = 0; channel < 4; ++channel) {
            auto& variable =
                problem_cfg.variables[
                    offset+state_size+2+channel];

            variable.free = true;

            variable.scale =
                (channel == 0 || channel == 2)
                    ? 1.0e6
                    : 10.0;
        }
    }

    // Keep constraint scaling explicit.
    const auto reference_constraints =
        assemble_full_space_trajectory_constraints(
            truth, cfg.parameters);

    problem_cfg.constraint_scales.assign(
        reference_constraints.combined().size(), 1.0);

    // All four boundary channels are observed without noise.
    // Sigma remains positive because it defines the weighting
    // of the least-squares objective.
    for (const auto& boundary : truth.reconstructed) {
        BoundaryMeasurement measurement;

        measurement.observed = boundary;
        measurement.available = {true, true, true, true};
        measurement.sigma = {
            1.0e3, 0.1, 1.0e3, 0.1
        };

        problem_cfg.measurements.push_back(measurement);
    }

    const FullSpaceProblem problem(
        std::move(problem_cfg));

    const auto zero =
        std::vector<double>(
            problem.free_variable_count(), 0.0);

    // The reference trajectory corresponds to u=0.
    const auto reference_evaluation =
        problem.evaluate(zero);

    require(
        reference_evaluation.measurement.value < 1.0e-20,
        "Reference measurement objective is nonzero");

    // Perturb the free optimization coordinates.
    // This leaves the fixed initial state and inputs intact.
    Eigen::VectorXd initial_coordinates(
        static_cast<Eigen::Index>(
            problem.free_variable_count()));

    for (Eigen::Index i = 0;
         i < initial_coordinates.size(); ++i) {

        initial_coordinates[i] =
            1.0e-4 *
            std::sin(0.7*static_cast<double>(i+1));
    }

    require(initial_coordinates.norm() > 0.0,
            "Initial perturbation vanished");

    DenseSQPOptions options;
    options.maximum_iterations = 40;
    options.feasibility_tolerance = 1.0e-8;
    options.stationarity_tolerance = 1.0e-8;

    const auto result =
        solve_dense_full_space_sqp(
            problem, initial_coordinates, options);

    std::cout
        << "SQP termination: "
        << static_cast<int>(result.termination) << '\n'
        << "Iterations recorded: "
        << result.history.size() << '\n'
        << "Final objective: "
        << result.final_objective << '\n'
        << "Final constraint L2: "
        << result.final_constraint_l2 << '\n'
        << "Final projected gradient: "
        << result.final_projected_gradient << '\n';

    for (const auto& iteration : result.history) {
        std::cout
            << "iter=" << iteration.iteration
            << " J=" << iteration.objective
            << " c2=" << iteration.constraint_l2
            << " pg=" << iteration.projected_gradient
            << " alpha=" << iteration.step_length
            << " rankA=" << iteration.constraint_rank
            << " rankB=" << iteration.reduced_measurement_rank
            << " unresolved=" << iteration.unresolved_directions
            << '\n';
    }

    require(result.converged(),
            "Synthetic pipeline reconstruction did not converge");

    std::vector<double> recovered_coordinates(
        result.coordinates.data(),
        result.coordinates.data() +
            result.coordinates.size());

    const auto reconstructed =
        problem.trajectory(recovered_coordinates);

    double maximum_pressure_error = 0.0;
    double maximum_flow_error = 0.0;

    for (std::size_t k = 0; k < levels; ++k) {
        for (std::size_t j = 0; j < state_size; ++j) {
            const double error = std::abs(
                reconstructed.states[k][j] -
                truth.states[k][j]);

            const bool pressure =
                j == 0 || (j >= 3 && j % 2 == 1);

            if (pressure)
                maximum_pressure_error =
                    std::max(maximum_pressure_error, error);
            else
                maximum_flow_error =
                    std::max(maximum_flow_error, error);
        }
    }

    const double reconstructed_defect =
        maximum_constraint_residual(
            reconstructed, cfg.parameters);

    std::cout
        << "Maximum pressure error [Pa]: "
        << maximum_pressure_error << '\n'
        << "Maximum flow error [kg/s]: "
        << maximum_flow_error << '\n'
        << "Reconstructed constraint Linf: "
        << reconstructed_defect << '\n'
        << "Forward mass defect [kg/s]: "
        << forward.summary.max_mass_defect << '\n';


    const auto conservation =
        audit_conservation(
            reconstructed, cfg.parameters);

    const auto rank =
        audit_rank(problem, result.coordinates);

    std::cout
        << "===== INDEPENDENT CONSERVATION AUDIT =====\n"
        << "Maximum physical mass imbalance [kg/s]: "
        << conservation.maximum_mass_imbalance << '\n'
        << "Maximum linepack change rate [kg/s]: "
        << conservation.maximum_linepack_change_rate << '\n'
        << "Maximum conservation defect [kg/s]: "
        << conservation.maximum_conservation_defect << '\n'
        << "Maximum integrated defect [kg]: "
        << conservation.maximum_integrated_defect << '\n';

    std::cout
        << "===== INDEPENDENT RANK AUDIT =====\n"
        << "Free variables: "
        << rank.variables << '\n'
        << "Equality constraints: "
        << rank.constraints << '\n'
        << "Measurement residuals: "
        << rank.measurements << '\n'
        << "Constraint rank: "
        << rank.constraint_rank << '\n'
        << "Constraint nullity: "
        << rank.nullity << '\n'
        << "Reduced measurement rank: "
        << rank.reduced_measurement_rank << '\n'
        << "Unresolved directions: "
        << rank.unresolved << '\n';

    require(rank.variables == 24,
            "Unexpected free-variable count");

    require(rank.constraints == 24,
            "Unexpected equality-constraint count");

    require(rank.constraint_rank == 24,
            "Constraint Jacobian is not full rank");

    require(rank.nullity == 0,
            "Unexpected constraint null space");

    require(rank.reduced_measurement_rank == 0,
            "Unexpected reduced measurement rank");

    require(rank.unresolved == 0,
            "Unexpected unresolved directions");

    require(
        conservation.maximum_conservation_defect < 1.0e-5,
        "Reconstructed mass-conservation defect too large");

    require(
        conservation.maximum_integrated_defect < 1.0e-5,
        "Integrated mass-conservation defect too large");

    require(result.final_constraint_l2 <=
                options.feasibility_tolerance,
            "Scaled nonlinear feasibility not achieved");

    // Stage A checks recovery of the known reference.
    // These are absolute physical-unit tolerances.
    require(maximum_pressure_error < 1.0,
            "Pressure reconstruction error too large");

    require(maximum_flow_error < 1.0e-4,
            "Flow reconstruction error too large");

    std::cout
        << "PASS E2.4b.11b.3 synthetic pipeline reconstruction\n";
}

} // namespace

int main()
{
    try {
        test_synthetic_reconstruction();
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E2.4b.11b.3: "
            << error.what() << '\n';

        return 1;
    }
}
