#include "tpsf/dense_full_space_sqp.hpp"
#include "tpsf/forecast_continuation.hpp"
#include "tpsf/dense_full_space_subproblem.hpp"
#include "tpsf/scheduled_forward.hpp"
#include "tpsf/pipeline_decision.hpp"
#include "tpsf/decision_timescales.hpp"
#include "tpsf/statefinding_quality.hpp"


#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
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
    const std::vector<std::size_t>& observation_levels,
    double pressure_sigma)
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
            std::find(
                observation_levels.begin(),
                observation_levels.end(),
                k) != observation_levels.end(),
            false
        };

        measurement.sigma = {
            1.0e3, 0.1, pressure_sigma, 0.1
        };

        cfg.measurements.push_back(measurement);
    }

    return FullSpaceProblem(std::move(cfg));
}


struct ExperimentCase {
    const char* name;
    std::vector<std::size_t> observations;
    double pressure_sigma;
};

void run_case(const ExperimentCase& experiment)
{
    const auto cfg = make_configuration();

    const InitialState initial{
        {5.0e6, 5.0e6, 5.0e6},
        {10.0, 10.0},
        10.0
    };

    const BoundarySchedule schedule{
        {5.0e6, 5.0e6, 5.0e6},
        {10.0, 10.1, 10.2}
    };

    const auto truth_forward =
        run_forward_scheduled_recorded(
            cfg, initial, schedule);

    require(
        truth_forward.summary.info == 0 &&
        truth_forward.summary.completed_steps == cfg.steps,
        "Truth forward integration failed");

    const auto truth =
        recorded_trajectory(cfg, truth_forward);

    InitialState guessed_initial = initial;

    guessed_initial.pressure.back() += 2.0e4;

    const auto guessed_forward =
        run_forward_scheduled_recorded(
            cfg, guessed_initial, schedule);

    require(
        guessed_forward.summary.info == 0 &&
        guessed_forward.summary.completed_steps == cfg.steps,
        "Independent reference integration failed");

    const auto reference =
        recorded_trajectory(cfg, guessed_forward);

    const auto problem =
        make_problem(
            reference,
            truth,
            cfg.parameters,
            experiment.observations,
            experiment.pressure_sigma);

    Eigen::VectorXd u0 = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(
            problem.free_variable_count()));

    DenseSQPOptions sqp_options;
    sqp_options.maximum_iterations = 40;
    sqp_options.feasibility_tolerance = 1.0e-8;
    sqp_options.stationarity_tolerance = 1.0e-8;

    const auto solution =
        solve_dense_full_space_sqp(
            problem, u0, sqp_options);

    // A missing-measurement solution can be feasible
    // without identifying the true initial condition.
    // Report that distinction rather than assuming
    // convergence implies reconstruction accuracy.

    const auto quality =
        assess_statefinding_quality(
            problem, solution);

    require(
        quality.feasible_dimension == 1,
        "Unexpected feasible dimension");

    require(
        quality.measurement_rank ==
            static_cast<std::size_t>(
                experiment.observations.empty() ? 0 : 1),
        "Unexpected measurement rank");

    require(
        quality.unresolved_directions ==
            static_cast<std::size_t>(
                experiment.observations.empty() ? 1 : 0),
        "Unexpected unresolved dimension");

    std::cout << std::scientific
              << std::setprecision(12);

    std::cout
        << "CASE " << experiment.name << '\n'
        << "  pressure_sigma_pa="
        << experiment.pressure_sigma << '\n'
        << "  sqp_converged=" << solution.converged() << '\n'
        << "  constraint_l2="
        << solution.final_constraint_l2 << '\n'
        << "  feasible_dimension="
        << quality.feasible_dimension << '\n'
        << "  measurement_rank="
        << quality.measurement_rank << '\n'
        << "  unresolved_dimension="
        << quality.unresolved_directions << '\n';

    if (quality.information_singular_values.size() > 0)
        std::cout
            << "  leading_information_singular_value="
            << quality.information_singular_values[0]
            << '\n';

    if (!solution.converged()) {
        std::cout
            << "  decision_status=not_evaluated\n";
        return;
    }

    const std::vector<double> recovered_coordinates(
        solution.coordinates.data(),
        solution.coordinates.data() +
            solution.coordinates.size());

    const auto reconstructed =
        problem.trajectory(recovered_coordinates);

    const double pressure_error =
        std::abs(
            reconstructed.states[0][5] -
            truth.states[0][5]);

    std::cout
        << "  initial_pressure_error_pa="
        << pressure_error << '\n';

    auto forecast_configuration = cfg;
    forecast_configuration.steps = 10;

    BoundarySchedule future;

    for (int k = 0;
         k <= forecast_configuration.steps;
         ++k) {

        const double time =
            k * forecast_configuration.parameters.dt;

        future.inlet_pressure.push_back(
            schedule.inlet_pressure.back());

        future.outlet_flow.push_back(
            schedule.outlet_flow.back() +
            0.02*time);
    }

    PipelineDecisionOptions decision_options;

    decision_options.minimum_outlet_pressure_pa =
        4.0e6;

    decision_options.measurement_radius = 1.0;

    // Do not introduce an artificial bound.
    decision_options.bound_unresolved = false;

    const auto decision =
        assess_pipeline_decision(
            problem,
            solution.coordinates,
            quality,
            forecast_configuration,
            future,
            decision_options);

    const auto timescales =
        assess_operational_timescales(decision);

    const auto& directional = decision.directional;

    // E5.7b: maximum admissible constant minimum-pressure
    // requirement over the sampled forecast.
    //
    // This is a local-linearized, sampled-time result.
    // It is not a nonlinear or continuous-time certificate.

    double critical_pressure_limit =
        std::numeric_limits<double>::infinity();

    bool finite_critical_limit = true;

    require(
        directional.time.size() ==
            directional.nominal_outlet_pressure.size(),
        "Forecast time/pressure dimension mismatch");

    require(
        directional.time.size() ==
            directional.measurement_unit_radius.size(),
        "Forecast measurement-radius dimension mismatch");

    require(
        directional.time.size() ==
            directional.unresolved_influence.size(),
        "Forecast unresolved-influence dimension mismatch");

    for (std::size_t k = 0;
         k < directional.time.size();
         ++k) {

        if (!decision_options.bound_unresolved &&
            directional.unresolved_influence[k]) {

            finite_critical_limit = false;
            break;
        }

        double radius =
            decision_options.measurement_radius *
            directional.measurement_unit_radius[k];

        if (decision_options.bound_unresolved) {
            radius +=
                decision_options.unresolved_radius *
                directional.unresolved_unit_radius[k];
        }

        const double lower_pressure =
            directional.nominal_outlet_pressure[k] -
            radius;

        critical_pressure_limit =
            std::min(
                critical_pressure_limit,
                lower_pressure);
    }

    if (finite_critical_limit) {
        std::cout
            << "  critical_pressure_limit_pa="
            << critical_pressure_limit << '\n';
    } else {
        std::cout
            << "  critical_pressure_limit_pa=UNBOUNDED_UNCERTAINTY\n";
    }


    // E5.7b: independent threshold-crossing validation.
    //
    // Recompute the engineering classification rather than
    // reusing the lower-envelope calculation.
    //
    // This validates the local-linearized decision machinery,
    // not the nonlinear forecast under finite perturbations.

    if (finite_critical_limit) {
        constexpr double pressure_offset = 1.0;

        auto below_options = decision_options;

        below_options.minimum_outlet_pressure_pa =
            critical_pressure_limit - pressure_offset;

        const auto below =
            assess_pipeline_decision(
                problem,
                solution.coordinates,
                quality,
                forecast_configuration,
                future,
                below_options);

        require(
            below.decision.classification ==
                DecisionClassification::certified_satisfactory,
            "Decision classifier rejected a subcritical limit");

        auto above_options = decision_options;

        above_options.minimum_outlet_pressure_pa =
            critical_pressure_limit + pressure_offset;

        const auto above =
            assess_pipeline_decision(
                problem,
                solution.coordinates,
                quality,
                forecast_configuration,
                future,
                above_options);

        require(
            above.decision.classification !=
                DecisionClassification::certified_satisfactory,
            "Decision classifier accepted a supercritical limit");

        std::cout
            << "  critical_boundary_validation=PASS\n";
    }
    else {
        // Unbounded unresolved uncertainty must not produce
        // a satisfactory forecast classification merely
        // because the engineering limit is lowered.

        auto low_options = decision_options;

        low_options.minimum_outlet_pressure_pa = 1.0e6;

        const auto low =
            assess_pipeline_decision(
                problem,
                solution.coordinates,
                quality,
                forecast_configuration,
                future,
                low_options);

        require(
            low.decision.classification !=
                DecisionClassification::certified_satisfactory,
            "Unbounded uncertainty produced a false certificate");

        std::cout
            << "  unbounded_uncertainty_validation=PASS\n";
    }

    const double maximum_measurement_radius =
        directional.measurement_unit_radius.empty()
            ? 0.0
            : *std::max_element(
                directional.measurement_unit_radius.begin(),
                directional.measurement_unit_radius.end());

    const double maximum_unresolved_radius =
        directional.unresolved_unit_radius.empty()
            ? 0.0
            : *std::max_element(
                directional.unresolved_unit_radius.begin(),
                directional.unresolved_unit_radius.end());

    std::cout
        << "  maximum_measurement_radius_pa="
        << maximum_measurement_radius << '\n'
        << "  maximum_unresolved_unit_radius_pa="
        << maximum_unresolved_radius << '\n'
        << "  satisfactory_samples="
        << timescales.decision.satisfactory_count << '\n'
        << "  violation_samples="
        << timescales.decision.violation_count << '\n'
        << "  indeterminate_samples="
        << timescales.decision.indeterminate_count << '\n';

    if (timescales.decision.resolved.time)
        std::cout
            << "  persistent_resolution_time_s="
            << *timescales.decision.resolved.time << '\n';
    else
        std::cout
            << "  persistent_resolution_time_s=NA\n";

    if (timescales.decision.satisfactory.time)
        std::cout
            << "  persistent_satisfactory_time_s="
            << *timescales.decision.satisfactory.time << '\n';
    else
        std::cout
            << "  persistent_satisfactory_time_s=NA\n";

    require(
        decision.attenuation.maximum_tangent_residual
            < 1.0e-6,
        "Excessive tangent residual");

    require(
        decision.attenuation.maximum_mass_defect
            < 1.0e-5,
        "Excessive mass conservation defect");

    std::cout
        << "PASS " << experiment.name << '\n';
}

} // namespace


int main()
{
    try {
        const std::vector<ExperimentCase> cases{
            {"none", {}, 1000.0},
            {"t0", {0}, 1000.0},
            {"t1", {1}, 1000.0},
            {"t2", {2}, 1000.0},
            {"all", {0, 1, 2}, 1000.0},
            {"t0_precise", {0}, 100.0},
            {"t1_precise", {1}, 100.0},
            {"t2_precise", {2}, 100.0},
            {"all_precise", {0, 1, 2}, 100.0}
        };

        for (const auto& experiment : cases)
            run_case(experiment);

        std::cout
            << "PASS E5.8a measurement design\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E5.7b: "
            << error.what() << '\n';

        return 1;
    }
}
