#include "tpsf/dense_full_space_sqp.hpp"
#include "tpsf/dense_full_space_subproblem.hpp"
#include "tpsf/forecast_continuation.hpp"
#include "tpsf/forecast_directional_uncertainty.hpp"
#include "tpsf/measurement_information.hpp"
#include "tpsf/pipeline_decision.hpp"
#include "tpsf/scheduled_forward.hpp"
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

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

Configuration configuration() {
    return {2, 2,
        pipe_sim::TransientParameters{100.0, 1.0, 0.65, 0.5,
            0.19634954084936207, 288.15, 0.9, 518.3,
            1.1e-5, 1.0e-5},
        pipe_sim::DemandRamp{5.0e6, 10.0, 10.2, 2.0},
        pipe_sim::NewtonOptions{1.0e-8, 1.0e-10, 30, false}};
}

FullSpaceTrajectory record(const Configuration& cfg,
                           const RecordedForwardResult& forward) {
    require(forward.summary.info == 0 &&
            forward.summary.completed_steps == cfg.steps,
            "Independent forward integration failed");
    FullSpaceTrajectory trajectory;
    trajectory.n = cfg.cells;
    for (const auto& snapshot : forward.trajectory.snapshots) {
        trajectory.states.push_back(snapshot.state);
        trajectory.inputs.push_back(snapshot.inputs);
        trajectory.reconstructed.push_back(snapshot.physical_channels);
    }
    require(trajectory.states.size() ==
            static_cast<std::size_t>(cfg.steps + 1),
            "Incomplete recorded trajectory");
    return trajectory;
}

FullSpaceProblem problem_for(const FullSpaceTrajectory& reference,
                             const FullSpaceTrajectory& truth,
                             const pipe_sim::TransientParameters& parameters,
                             const std::vector<std::size_t>& observed) {
    FullSpaceProblemConfiguration cfg;
    cfg.reference = reference;
    cfg.parameters = parameters;
    const std::size_t state_size = static_cast<std::size_t>(2*reference.n+2);
    const std::size_t block = state_size + 6;
    const std::size_t levels = reference.states.size();
    cfg.variables.assign(levels*block, ProblemVariable{});
    for (std::size_t k=0; k<levels; ++k) {
        const auto offset = k*block;
        for (std::size_t j=0; j<state_size; ++j) {
            auto& variable = cfg.variables[offset+j];
            variable.free = (k != 0 || j == 3 || j == 5);
            variable.scale = (j == 0 || (j >= 3 && j%2 == 1))
                ? 1.0e6 : 10.0;
        }
        cfg.variables[offset+state_size].free = false;
        cfg.variables[offset+state_size+1].free = false;
        for (std::size_t channel=0; channel<4; ++channel) {
            auto& variable = cfg.variables[offset+state_size+2+channel];
            variable.free = true;
            variable.scale = (channel == 0 || channel == 2)
                ? 1.0e6 : 10.0;
        }
    }
    const auto constraints =
        assemble_full_space_trajectory_constraints(reference, parameters);
    cfg.constraint_scales.assign(constraints.combined().size(), 1.0);
    for (std::size_t k=0; k<levels; ++k) {
        BoundaryMeasurement measurement;
        measurement.observed = truth.reconstructed[k];
        measurement.available = {false, false,
            std::find(observed.begin(), observed.end(), k) != observed.end(),
            false};
        measurement.sigma = {1.0e3, 0.1, 1.0e3, 0.1};
        cfg.measurements.push_back(measurement);
    }
    return FullSpaceProblem(std::move(cfg));
}

void run_case(const char* label,
              const std::vector<std::size_t>& observed,
              const FullSpaceTrajectory& truth,
              const FullSpaceTrajectory& reference,
              const Configuration& cfg,
              const BoundarySchedule& schedule) {
    auto problem = problem_for(reference, truth, cfg.parameters, observed);
    Eigen::VectorXd initial = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(problem.free_variable_count()));
    DenseSQPOptions options;
    options.maximum_iterations = 40;
    options.feasibility_tolerance = 1.0e-8;
    options.stationarity_tolerance = 1.0e-8;
    const auto solution = solve_dense_full_space_sqp(problem, initial, options);
    const auto quality = assess_statefinding_quality(problem, solution);
    const std::size_t expected_rank = observed.size();
    std::cout << "CASE " << label << '\n'
              << "  converged=" << solution.converged() << '\n'
              << "  constraint_l2=" << solution.final_constraint_l2 << '\n'
              << "  feasible_dimension=" << quality.feasible_dimension << '\n'
              << "  measurement_rank=" << quality.measurement_rank << '\n'
              << "  unresolved_directions=" << quality.unresolved_directions << '\n'
              << "  information_singular_values="
              << quality.information_singular_values.transpose() << '\n';
    require(quality.feasible_dimension == 2,
            "Expected two feasible initial-pressure directions");
    require(quality.measurement_rank == expected_rank,
            "Unexpected reduced measurement rank");
    require(quality.unresolved_directions == 2-expected_rank,
            "Unexpected unresolved dimension");
    require(solution.converged(), "SQP reconstruction did not converge");

    const std::vector<double> coordinates(
        solution.coordinates.data(),
        solution.coordinates.data()+solution.coordinates.size());
    const auto reconstructed = problem.trajectory(coordinates);
    std::cout << "  interior_initial_error_pa="
              << reconstructed.states[0][3]-truth.states[0][3] << '\n'
              << "  outlet_initial_error_pa="
              << reconstructed.states[0][5]-truth.states[0][5] << '\n';

    auto forecast_cfg = cfg;
    forecast_cfg.steps = 10;
    BoundarySchedule future;
    for (int k=0; k<=forecast_cfg.steps; ++k) {
        const double t = k*forecast_cfg.parameters.dt;
        future.inlet_pressure.push_back(schedule.inlet_pressure.back());
        future.outlet_flow.push_back(schedule.outlet_flow.back()+0.02*t);
    }
    PipelineDecisionOptions decision_options;
    decision_options.minimum_outlet_pressure_pa = 4.0e6;
    decision_options.measurement_radius = 1.0;
    decision_options.bound_unresolved = false;
    const auto decision = assess_pipeline_decision(
        problem, solution.coordinates, quality, forecast_cfg,
        future, decision_options);
    const auto& d = decision.directional;
    // E5.9a: independent nonlinear finite-amplitude validation.
    // The directions are scaled full-space coordinates. Map their
    // terminal physical-state increments into the original solver.
    // Do not modify future boundary controls.
    if (expected_rank == 1) {
        const auto assembled = assemble_dense_full_space_subproblem(
            problem, solution.coordinates);
        const auto information = decompose_measurement_information(
            assembled.measurement_jacobian, quality.feasible_basis,
            quality.relative_rank_tolerance);
        require(information.resolved_basis.cols() == 1 &&
                information.unresolved_basis.cols() == 1,
                "Expected one resolved and one unresolved direction");
        const auto nominal_initial = initial_state_from_full_space(
            reconstructed, reconstructed.states.size()-1);
        const auto nominal_forward = run_forward_scheduled(
            forecast_cfg, nominal_initial, future);
        require(nominal_forward.info == 0,
                "Nonlinear nominal forecast failed");

        auto check_direction = [&](const char* name,
                                   const Eigen::VectorXd& basis,
                                   const Eigen::VectorXd& predicted) {
            const std::vector<double> basis_coordinates(
                basis.data(), basis.data()+basis.size());
            const auto mapped = problem.trajectory_direction(
                basis_coordinates);
            const auto& physical = mapped.states.back();
            require(physical.size() == 6,
                    "Unexpected terminal physical state size");
            require(std::abs(physical[0]) < 1.0e-7,
                    "Direction changes fixed inlet pressure");
            require(predicted.size() ==
                    static_cast<Eigen::Index>(nominal_forward.outlet_pressure.size()),
                    "Nonlinear forecast output length mismatch");

            auto perturb = [&](double h) {
                InitialState state = nominal_initial;
                state.pressure[0] += h*physical[0];
                state.inlet_flow += h*physical[1];
                for (int i=0; i<2; ++i) {
                    state.flow[i] += h*physical[2*i+2];
                    state.pressure[i+1] += h*physical[2*i+3];
                }
                return state;
            };

            for (double h : {1.0e-2, 1.0e-3, 1.0e-4}) {
                const auto plus = run_forward_scheduled(
                    forecast_cfg, perturb(h), future);
                const auto minus = run_forward_scheduled(
                    forecast_cfg, perturb(-h), future);
                require(plus.info == 0 && minus.info == 0,
                        "Perturbed nonlinear forecast failed");
                double derivative_error = 0.0;
                double remainder = 0.0;
                double sensitivity_scale = 0.0;
                for (Eigen::Index k=0; k<predicted.size(); ++k) {
                    const auto j = static_cast<std::size_t>(k);
                    const double derivative =
                        (plus.outlet_pressure[j]-minus.outlet_pressure[j])/(2*h);
                    derivative_error = std::max(derivative_error,
                        std::abs(derivative-predicted[k]));
                    remainder = std::max(remainder, std::abs(
                        plus.outlet_pressure[j]-nominal_forward.outlet_pressure[j]
                        -h*predicted[k]));
                    sensitivity_scale = std::max(sensitivity_scale,
                        std::abs(predicted[k]));
                }
                // Small-h cancellation limits absolute derivative accuracy.
                // The tolerance is deliberately reported alongside raw errors.
                require(derivative_error <=
                        2.0e-2+1.0e-3*sensitivity_scale,
                        "Nonlinear centered derivative disagrees with tangent");
                std::cout << "  e59a_" << name << "_h=" << h
                          << " derivative_error_pa_per_coordinate="
                          << derivative_error
                          << " finite_remainder_pa=" << remainder << '\n';
            }
        };
        check_direction("resolved", information.resolved_basis.col(0),
                        d.resolved_sensitivity.col(0));
        check_direction("unresolved", information.unresolved_basis.col(0),
                        d.unresolved_sensitivity.col(0));
        std::cout << "  e59a_nonlinear_directional_validation=PASS\n";
    }
    require(d.unresolved_sensitivity.cols() ==
            static_cast<Eigen::Index>(quality.unresolved_directions),
            "Unresolved sensitivity dimension mismatch");
    for (Eigen::Index j=0; j<d.unresolved_sensitivity.cols(); ++j) {
        double maximum=0.0;
        for (Eigen::Index k=0; k<d.unresolved_sensitivity.rows(); ++k)
            maximum=std::max(maximum,std::abs(d.unresolved_sensitivity(k,j)));
        std::cout << "  unresolved_direction_" << j
                  << "_max_sensitivity_pa=" << maximum << '\n';
    }
    // E5.7c: explicitly bounded unresolved uncertainty.
    // Radius is Euclidean displacement in scaled feasible coordinates,
    // not a statistical confidence interval or physical pressure bound.
    if (expected_rank == 1) {
        require(d.time.size() == d.nominal_outlet_pressure.size() &&
                d.time.size() == d.measurement_unit_radius.size() &&
                d.time.size() == d.unresolved_unit_radius.size(),
                "Directional forecast size mismatch");

        const double limit = decision_options.minimum_outlet_pressure_pa;
        double critical_radius = std::numeric_limits<double>::infinity();
        double minimum_margin = std::numeric_limits<double>::infinity();
        bool has_sensitive_sample = false;

        for (std::size_t k=0; k<d.time.size(); ++k) {
            const double margin = d.nominal_outlet_pressure[k]
                - decision_options.measurement_radius
                    * d.measurement_unit_radius[k] - limit;
            const double sensitivity = d.unresolved_unit_radius[k];
            require(std::isfinite(margin) && std::isfinite(sensitivity)
                    && sensitivity >= 0.0,
                    "Invalid radius boundary inputs");
            minimum_margin = std::min(minimum_margin, margin);
            if (sensitivity > 0.0) {
                has_sensitive_sample = true;
                critical_radius = std::min(critical_radius,
                                           margin / sensitivity);
            } else {
                require(margin >= 0.0,
                        "Zero-sensitivity sample violates pressure limit");
            }
        }
        require(has_sensitive_sample && minimum_margin > 0.0 &&
                std::isfinite(critical_radius) && critical_radius > 0.0,
                "No positive bounded-uncertainty decision boundary");

        // Independent classifier calls on either side of the boundary.
        // The boundary is computed from the directional forecast arrays,
        // not from the classifier's reported aggregate classification.
        const double offset = std::max(1.0e-3, 1.0e-3*critical_radius);
        require(offset < critical_radius,
                "Boundary perturbation exceeds critical radius");
        auto bounded = decision_options;
        bounded.bound_unresolved = true;
        bounded.unresolved_radius = critical_radius - offset;
        const auto below = assess_pipeline_decision(
            problem, solution.coordinates, quality, forecast_cfg,
            future, bounded);
        require(below.decision.classification ==
                    DecisionClassification::certified_satisfactory,
                "Subcritical bounded uncertainty was not satisfactory");

        bounded.unresolved_radius = critical_radius + offset;
        const auto above = assess_pipeline_decision(
            problem, solution.coordinates, quality, forecast_cfg,
            future, bounded);
        require(above.decision.classification !=
                    DecisionClassification::certified_satisfactory,
                "Supercritical bounded uncertainty was satisfactory");

        bounded.unresolved_radius = 1.0;
        const auto unit = assess_pipeline_decision(
            problem, solution.coordinates, quality, forecast_cfg,
            future, bounded);
        require(unit.decision.classification ==
                    DecisionClassification::certified_satisfactory,
                "Unit bounded uncertainty was not satisfactory");

        require(decision.decision.classification !=
                    DecisionClassification::certified_satisfactory,
                "Unbounded unresolved uncertainty was certified satisfactory");

        std::cout << "  pressure_limit_pa=" << limit << '\n'
                  << "  minimum_zero_unresolved_margin_pa="
                  << minimum_margin << '\n'
                  << "  critical_unresolved_coordinate_radius="
                  << critical_radius << '\n'
                  << "  subcritical_radius=" << critical_radius-offset << '\n'
                  << "  supercritical_radius=" << critical_radius+offset << '\n'
                  << "  unit_bounded_classification="
                  << static_cast<int>(unit.decision.classification) << '\n'
                  << "  bounded_boundary_validation=PASS\n";
    }
    std::cout << "  unbounded_decision_classification="
              << static_cast<int>(decision.decision.classification) << '\n'
              << "  tangent_residual=" << d.maximum_tangent_residual << '\n'
              << "  tangent_mass_defect=" << d.maximum_tangent_mass_defect << '\n'
              << "  nominal_mass_defect=" << d.maximum_nominal_mass_defect << '\n';
}
}

int main() {
    try {
        std::cout << std::scientific << std::setprecision(12);
        const auto cfg = configuration();
        const InitialState truth_initial{{5.0e6,5.0e6,5.0e6},
                                         {10.0,10.0},10.0};
        const BoundarySchedule schedule{{5.0e6,5.0e6,5.0e6},
                                        {10.0,10.1,10.2}};
        const auto truth = record(cfg,
            run_forward_scheduled_recorded(cfg,truth_initial,schedule));
        InitialState guessed = truth_initial;
        guessed.pressure[1] += 1.0e4;
        guessed.pressure.back() += 2.0e4;
        const auto reference = record(cfg,
            run_forward_scheduled_recorded(cfg,guessed,schedule));
        run_case("none",{},truth,reference,cfg,schedule);
        run_case("t0",{0},truth,reference,cfg,schedule);
        run_case("t0_t1",{0,1},truth,reference,cfg,schedule);
        std::cout << "PASS E5.7c structural experiment\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL E5.7c: " << e.what() << '\n';
        return 1;
    }
}
