#include "tpsf/dense_full_space_sqp.hpp"
#include "tpsf/forecast_continuation.hpp"
#include "tpsf/dense_full_space_subproblem.hpp"
#include "tpsf/scheduled_forward.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace tpsf;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

Configuration make_configuration()
{
    return {
        2,
        2,
        pipe_sim::TransientParameters{
            100.0,
            1.0,
            0.65,
            0.5,
            0.19634954084936207,
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

FullSpaceTrajectory recorded_trajectory(
    const Configuration& cfg,
    const RecordedForwardResult& forward)
{
    FullSpaceTrajectory trajectory;
    trajectory.n = cfg.cells;

    for (const auto& snapshot : forward.trajectory.snapshots) {
        trajectory.states.push_back(snapshot.state);
        trajectory.inputs.push_back(snapshot.inputs);
        trajectory.reconstructed.push_back(
            snapshot.physical_channels);
    }

    require(
        trajectory.states.size() ==
            static_cast<std::size_t>(cfg.steps+1),
        "Incomplete forward trajectory");

    return trajectory;
}

Eigen::Index numerical_rank(
    const Eigen::VectorXd& singular_values)
{
    if (singular_values.size() == 0 ||
        singular_values[0] == 0.0)
        return 0;

    constexpr double relative_tolerance = 1.0e-12;

    Eigen::Index rank = 0;

    for (Eigen::Index i = 0;
         i < singular_values.size(); ++i) {

        if (singular_values[i] >
            relative_tolerance*singular_values[0])
            ++rank;
    }

    return rank;
}

struct RankAudit {
    Eigen::Index variables = 0;
    Eigen::Index constraints = 0;
    Eigen::Index measurements = 0;
    Eigen::Index constraint_rank = 0;
    Eigen::Index nullity = 0;
    Eigen::Index measurement_rank = 0;
    Eigen::Index unresolved = 0;
    double information_singular_value = 0.0;
};

RankAudit audit_rank(
    const FullSpaceProblem& problem,
    const Eigen::VectorXd& coordinates)
{
    const auto assembled =
        assemble_dense_full_space_subproblem(
            problem, coordinates);

    const auto& A = assembled.constraint_jacobian;
    const auto& B = assembled.measurement_jacobian;

    RankAudit audit;

    audit.variables = A.cols();
    audit.constraints = A.rows();
    audit.measurements = B.rows();

    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        A, Eigen::ComputeFullV);

    audit.constraint_rank =
        numerical_rank(svd.singularValues());

    audit.nullity = A.cols()-audit.constraint_rank;

    if (audit.nullity == 0)
        return audit;

    const Eigen::MatrixXd Z =
        svd.matrixV().rightCols(audit.nullity);

    const Eigen::MatrixXd reduced = B*Z;

    // An empty measurement matrix has rank zero.
    // Eigen::JacobiSVD requires nonempty dimensions.
    if (reduced.rows() == 0 || reduced.cols() == 0) {
        audit.measurement_rank = 0;
        audit.unresolved = audit.nullity;
        return audit;
    }

    const Eigen::JacobiSVD<Eigen::MatrixXd> reduced_svd(
        reduced, Eigen::ComputeThinU |
                 Eigen::ComputeThinV);

    audit.information_singular_value =
        reduced_svd.singularValues()[0];

    audit.measurement_rank =
        numerical_rank(reduced_svd.singularValues());

    audit.unresolved =
        audit.nullity-audit.measurement_rank;

    return audit;
}

FullSpaceProblem make_problem(
    const FullSpaceTrajectory& reference,
    const FullSpaceTrajectory& truth,
    const pipe_sim::TransientParameters& parameters,
    bool observe_outlet_pressure,
    std::size_t observation_level)
{
    FullSpaceProblemConfiguration cfg;

    cfg.reference = reference;
    cfg.parameters = parameters;

    const std::size_t levels = reference.states.size();
    const std::size_t state_size =
        static_cast<std::size_t>(2*reference.n+2);

    const std::size_t block_size = state_size+6;

    cfg.variables.assign(
        levels*block_size, ProblemVariable{});

    for (std::size_t k = 0; k < levels; ++k) {
        const std::size_t offset = k*block_size;

        for (std::size_t j = 0; j < state_size; ++j) {
            auto& variable = cfg.variables[offset+j];

            variable.free = (k != 0 || j == 5);

            const bool pressure =
                j == 0 || (j >= 3 && j % 2 == 1);

            variable.scale =
                pressure ? 1.0e6 : 10.0;
        }

        // Prescribed boundary controls are fixed.
        cfg.variables[offset+state_size].free = false;
        cfg.variables[offset+state_size+1].free = false;

        // Reconstructed boundary channels remain free.
        for (std::size_t channel = 0;
             channel < 4; ++channel) {

            auto& variable =
                cfg.variables[
                    offset+state_size+2+channel];

            variable.free = true;

            variable.scale =
                (channel == 0 || channel == 2)
                    ? 1.0e6
                    : 10.0;
        }
    }

    const auto constraints =
        assemble_full_space_trajectory_constraints(
            reference, parameters);

    cfg.constraint_scales.assign(
        constraints.combined().size(), 1.0);

    for (std::size_t k = 0; k < levels; ++k) {
        BoundaryMeasurement measurement;

        measurement.observed =
            truth.reconstructed[k];

        measurement.available = {
            false,
            false,
            observe_outlet_pressure &&
                k == observation_level,
            false
        };

        measurement.sigma = {
            1.0e3, 0.1, 1.0e3, 0.1
        };

        cfg.measurements.push_back(measurement);
    }

    return FullSpaceProblem(std::move(cfg));
}

void test_independent_initial_pressure(std::size_t observation_level)
{
    const auto cfg = make_configuration();

    require(
        observation_level <= static_cast<std::size_t>(cfg.steps),
        "Observation level exceeds simulation horizon");

    std::cout
        << "Observation level: " << observation_level
        << ", time [s]: "
        << observation_level * cfg.parameters.dt
        << '\n';

    const InitialState initial{
        {5.0e6, 5.0e6, 5.0e6},
        {10.0, 10.0},
        10.0
    };

    const BoundarySchedule schedule{
        {5.0e6, 5.0e6, 5.0e6},
        {10.0, 10.1, 10.2}
    };

    const auto forward =
        run_forward_scheduled_recorded(
            cfg, initial, schedule);

    require(forward.summary.info == 0,
            "Forward solver failed");

    require(forward.summary.completed_steps == cfg.steps,
            "Incomplete forward simulation");

    const auto truth =
        recorded_trajectory(cfg, forward);

    // Deliberately incorrect initial outlet pressure.
    // This is a local-identification test; the remaining
    // trajectory reference is still forward-generated truth.
    constexpr double initial_pressure_error = 2.0e4;

    InitialState guessed_initial = initial;

    guessed_initial.pressure.back() +=
        initial_pressure_error;

    const auto guessed_forward =
        run_forward_scheduled_recorded(
            cfg, guessed_initial, schedule);

    require(
        guessed_forward.summary.info == 0,
        "Independent guessed forward solve failed");

    require(
        guessed_forward.summary.completed_steps ==
            cfg.steps,
        "Independent guessed forward solve incomplete");

    const auto reference =
        recorded_trajectory(cfg, guessed_forward);

    require(
        std::abs(
            reference.states[0][5] -
            truth.states[0][5] -
            initial_pressure_error) < 1.0e-8,
        "Independent initial pressure offset incorrect");

    const double reference_constraint_error =
        [&]() {
            const auto residuals =
                assemble_full_space_trajectory_constraints(
                    reference, cfg.parameters);

            double maximum = 0.0;

            for (double value : residuals.combined())
                maximum =
                    std::max(maximum, std::abs(value));

            return maximum;
        }();

    std::cout
        << "Independent reference constraint Linf: "
        << reference_constraint_error << '\n';

    require(
        reference_constraint_error < 1.0e-5,
        "Independent forward reference violates constraints");

    const auto problem =
        make_problem(
            reference,
            truth,
            cfg.parameters,
            true,
            observation_level);

    const auto unobserved_problem =
        make_problem(
            reference,
            truth,
            cfg.parameters,
            false,
            observation_level);

    Eigen::VectorXd u0 =
        Eigen::VectorXd::Zero(
            static_cast<Eigen::Index>(
                problem.free_variable_count()));

    const auto initial_rank =
        audit_rank(problem, u0);

    const auto missing_rank =
        audit_rank(unobserved_problem, u0);

    std::cout << std::scientific
              << std::setprecision(12);

    std::cout
        << "===== INITIAL RANK AUDIT =====\n"
        << "Free variables: "
        << initial_rank.variables << '\n'
        << "Equality constraints: "
        << initial_rank.constraints << '\n'
        << "Measurement residuals: "
        << initial_rank.measurements << '\n'
        << "Constraint rank: "
        << initial_rank.constraint_rank << '\n'
        << "Constraint nullity: "
        << initial_rank.nullity << '\n'
        << "Reduced measurement rank: "
        << initial_rank.measurement_rank << '\n'
        << "Unresolved directions: "
        << initial_rank.unresolved << '\n';

    std::cout
        << "===== MISSING-MEASUREMENT CONTROL =====\n"
        << "Reduced measurement rank: "
        << missing_rank.measurement_rank << '\n'
        << "Unresolved directions: "
        << missing_rank.unresolved << '\n';

    require(initial_rank.variables == 25,
            "Expected 25 free variables");

    require(initial_rank.constraints == 24,
            "Expected 24 equality constraints");

    require(initial_rank.constraint_rank == 24,
            "Unexpected constraint rank");

    require(initial_rank.nullity == 1,
            "Expected one feasible direction");

    require(initial_rank.measurements == 1,
            "Expected one measurement residual");

    require(initial_rank.measurement_rank == 1,
            "Outlet pressure is not locally identifiable");

    require(initial_rank.unresolved == 0,
            "Unexpected unresolved direction");

    require(missing_rank.nullity == 1,
            "Missing-measurement nullity changed");

    require(missing_rank.measurement_rank == 0,
            "Unexpected information without measurements");

    require(missing_rank.unresolved == 1,
            "Missing measurement should leave one direction");

    DenseSQPOptions options;

    options.maximum_iterations = 40;
    options.feasibility_tolerance = 1.0e-8;
    options.stationarity_tolerance = 1.0e-8;

    const auto result =
        solve_dense_full_space_sqp(
            problem, u0, options);

    std::cout
        << "===== RECONSTRUCTION =====\n"
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
            << " unresolved="
            << iteration.unresolved_directions
            << '\n';
    }

    require(result.converged(),
            "Unknown initial pressure reconstruction failed");

    const std::vector<double> recovered_coordinates(
        result.coordinates.data(),
        result.coordinates.data() +
            result.coordinates.size());

    const auto reconstructed =
        problem.trajectory(recovered_coordinates);

    const double recovered_pressure =
        reconstructed.states[0][5];

    const double true_pressure =
        truth.states[0][5];

    const double pressure_error =
        std::abs(recovered_pressure-true_pressure);

    const auto final_rank =
        audit_rank(problem, result.coordinates);

    std::cout
        << "True initial outlet pressure [Pa]: "
        << true_pressure << '\n'
        << "Initial pressure guess [Pa]: "
        << reference.states[0][5] << '\n'
        << "Recovered initial outlet pressure [Pa]: "
        << recovered_pressure << '\n'
        << "Recovered pressure error [Pa]: "
        << pressure_error << '\n'
        << "Final constraint rank: "
        << final_rank.constraint_rank << '\n'
        << "Final reduced measurement rank: "
        << final_rank.measurement_rank << '\n';

    require(
        result.final_constraint_l2 <=
            options.feasibility_tolerance,
        "Nonlinear feasibility not achieved");

    require(pressure_error < 1.0,
            "Initial outlet pressure recovery error too large");

    require(final_rank.constraint_rank == 24,
            "Final constraint rank changed");

    require(final_rank.measurement_rank == 1,
            "Final measurement rank changed");

    std::cout
        << "E5.7 initial information singular value: "
        << initial_rank.information_singular_value << '\n'
        << "E5.7 final information singular value: "
        << final_rank.information_singular_value << '\n';


    // --------------------------------------------------
    // E3: Continue the independently reconstructed state.
    // --------------------------------------------------

    auto forecast_configuration = cfg;
    forecast_configuration.steps = 10;

    BoundarySchedule future;

    const double terminal_time =
        static_cast<double>(cfg.steps) *
        cfg.parameters.dt;

    const double terminal_outlet_flow =
        schedule.outlet_flow.back();

    for (int k = 0; k <= forecast_configuration.steps; ++k) {
        const double future_time =
            static_cast<double>(k) *
            forecast_configuration.parameters.dt;

        future.inlet_pressure.push_back(
            schedule.inlet_pressure.back());

        // Continuous at the reconstruction/forecast
        // junction, followed by a modest demand increase.
        future.outlet_flow.push_back(
            terminal_outlet_flow +
            0.02*future_time);
    }

    const auto forecast =
        forecast_from_full_space(
            forecast_configuration,
            reconstructed,
            future);

    const auto reference_forecast =
        forecast_from_full_space(
            forecast_configuration,
            truth,
            future);

    require(
        forecast.summary.info == 0 &&
        forecast.summary.completed_steps ==
            forecast_configuration.steps,
        "SQP-initialized forecast failed");

    require(
        reference_forecast.summary.info == 0 &&
        reference_forecast.summary.completed_steps ==
            forecast_configuration.steps,
        "Reference forecast failed");

    double maximum_pressure_error = 0.0;
    double maximum_linepack_error = 0.0;

    for (int k = 0; k <= forecast_configuration.steps; ++k) {
        const auto index = static_cast<std::size_t>(k);

        maximum_pressure_error = std::max(
            maximum_pressure_error,
            std::abs(
                forecast.summary.outlet_pressure[index] -
                reference_forecast.summary.outlet_pressure[index]));

        maximum_linepack_error = std::max(
            maximum_linepack_error,
            std::abs(
                forecast.summary.linepack[index] -
                reference_forecast.summary.linepack[index]));
    }

    std::cout
        << "===== SQP-INITIALIZED FORECAST =====\n"
        << "Reconstruction terminal time [s]: "
        << terminal_time << '\n'
        << "Forecast duration [s]: "
        << forecast_configuration.steps *
           forecast_configuration.parameters.dt << '\n'
        << "Maximum outlet pressure error [Pa]: "
        << maximum_pressure_error << '\n'
        << "Maximum linepack error [kg]: "
        << maximum_linepack_error << '\n'
        << "Maximum conservation defect [kg/s]: "
        << forecast.summary.max_mass_defect << '\n';

    require(
        maximum_pressure_error < 1.0,
        "SQP forecast outlet-pressure error exceeds 1 Pa");

    require(
        maximum_linepack_error < 1.0e-3,
        "SQP forecast linepack error exceeds 1 g");

    require(
        forecast.summary.max_mass_defect < 1.0e-5,
        "SQP forecast mass-conservation defect too large");

    std::cout
        << "PASS E3 SQP reconstruction-to-forecast integration\n";

    std::cout
        << "PASS E2.4b.11b.4b independent initial pressure\n";
}

} // namespace


int main(int argc, char** argv)
{
    try {
        require(argc == 2,
                "Usage: independent_initial_pressure_test "
                "<observation-level>");

        const std::string argument = argv[1];

        require(argument == "0" ||
                argument == "1" ||
                argument == "2",
                "Observation level must be 0, 1, or 2");

        const std::size_t observation_level =
            static_cast<std::size_t>(std::stoul(argument));

        test_independent_initial_pressure(observation_level);

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E2.4b.11b.4b: "
            << error.what() << '\n';

        return 1;
    }
}
