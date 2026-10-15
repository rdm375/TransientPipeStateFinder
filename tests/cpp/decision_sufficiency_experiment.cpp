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

        // E5.9b: finite-amplitude nonlinear validation along the
        // unresolved direction. Future boundary inputs remain fixed.
        {
            const Eigen::VectorXd basis =
                information.unresolved_basis.col(0);
            const std::vector<double> basis_coordinates(
                basis.data(), basis.data() + basis.size());

            const auto mapped =
                problem.trajectory_direction(basis_coordinates);
            const auto& physical = mapped.states.back();

            const Eigen::VectorXd predicted =
                d.unresolved_sensitivity.col(0);

            require(physical.size() == 6,
                    "Unexpected E5.9b physical direction size");

            const double pressure_limit =
                decision_options.minimum_outlet_pressure_pa;

            const double measurement_radius =
                decision_options.measurement_radius;

            double critical_radius =
                std::numeric_limits<double>::infinity();

            for (std::size_t k = 0; k < d.time.size(); ++k) {
                const double margin =
                    d.nominal_outlet_pressure[k]
                    - measurement_radius * d.measurement_unit_radius[k]
                    - pressure_limit;

                const double sensitivity =
                    d.unresolved_unit_radius[k];

                if (sensitivity > 0.0)
                    critical_radius =
                        std::min(critical_radius, margin / sensitivity);
            }

            require(std::isfinite(critical_radius) &&
                    critical_radius > 0.0,
                    "Invalid E5.9b critical radius");

            std::cout
                << "  e59b_linearized_critical_radius="
                << critical_radius << '\n';

            auto perturb = [&](double h) {
                InitialState state = nominal_initial;

                state.pressure[0] += h * physical[0];
                state.inlet_flow += h * physical[1];

                for (int i = 0; i < 2; ++i) {
                    state.flow[i] += h * physical[2*i + 2];
                    state.pressure[i+1] += h * physical[2*i + 3];
                }

                return state;
            };


            // E5.9c: inspect the full pressure history at the
            // linearized critical unresolved-coordinate radius.
            //
            // This is a directional nonlinear diagnostic, not a
            // certificate over the complete uncertainty set.
            {
                std::cout
                    << "  e59c_active_constraint_diagnostic=BEGIN\n";

                for (double sign : {-1.0, 1.0}) {
                    const double h = sign * critical_radius;

                    const auto nonlinear =
                        run_forward_scheduled(
                            forecast_cfg, perturb(h), future);

                    require(nonlinear.info == 0,
                            "E5.9c nonlinear forecast failed");

                    require(nonlinear.outlet_pressure.size() ==
                            d.time.size(),
                            "E5.9c forecast length mismatch");

                    double minimum_linear_margin =
                        std::numeric_limits<double>::infinity();

                    double minimum_nonlinear_margin =
                        std::numeric_limits<double>::infinity();

                    std::size_t active_linear_index = 0;
                    std::size_t active_nonlinear_index = 0;
                    std::size_t maximum_remainder_index = 0;

                    double maximum_absolute_remainder = 0.0;

                    for (std::size_t k = 0;
                         k < d.time.size(); ++k) {

                        const double nominal =
                            nominal_forward.outlet_pressure[k];

                        const double linear =
                            nominal + h * predicted[
                                static_cast<Eigen::Index>(k)];

                        const double actual =
                            nonlinear.outlet_pressure[k];

                        const double remainder =
                            actual - linear;

                        const double allowance =
                            measurement_radius
                            * d.measurement_unit_radius[k];

                        const double linear_margin =
                            linear - allowance - pressure_limit;

                        const double nonlinear_margin =
                            actual - allowance - pressure_limit;

                        if (linear_margin < minimum_linear_margin) {
                            minimum_linear_margin = linear_margin;
                            active_linear_index = k;
                        }

                        if (nonlinear_margin <
                            minimum_nonlinear_margin) {
                            minimum_nonlinear_margin =
                                nonlinear_margin;
                            active_nonlinear_index = k;
                        }

                        if (std::abs(remainder) >
                            maximum_absolute_remainder) {
                            maximum_absolute_remainder =
                                std::abs(remainder);
                            maximum_remainder_index = k;
                        }

                        std::cout
                            << "  e59c_h=" << h
                            << " k=" << k
                            << " time=" << d.time[k]
                            << " nominal_pa=" << nominal
                            << " linear_pa=" << linear
                            << " nonlinear_pa=" << actual
                            << " remainder_pa=" << remainder
                            << " measurement_allowance_pa="
                            << allowance
                            << " linear_margin_pa="
                            << linear_margin
                            << " nonlinear_margin_pa="
                            << nonlinear_margin
                            << '\n';
                    }

                    const std::size_t a = active_linear_index;

                    const double active_linear_pressure =
                        nominal_forward.outlet_pressure[a]
                        + h * predicted[
                            static_cast<Eigen::Index>(a)];

                    const double active_nonlinear_pressure =
                        nonlinear.outlet_pressure[a];

                    std::cout
                        << "  e59c_summary_h=" << h
                        << " active_linear_index="
                        << active_linear_index
                        << " active_nonlinear_index="
                        << active_nonlinear_index
                        << " maximum_remainder_index="
                        << maximum_remainder_index
                        << " active_remainder_pa="
                        << active_nonlinear_pressure
                           - active_linear_pressure
                        << " maximum_absolute_remainder_pa="
                        << maximum_absolute_remainder
                        << " minimum_linear_margin_pa="
                        << minimum_linear_margin
                        << " minimum_nonlinear_margin_pa="
                        << minimum_nonlinear_margin
                        << '\n';
                }

                std::cout
                    << "  e59c_active_constraint_diagnostic=COMPLETE\n";
            }

            const std::vector<double> amplitudes = {
                0.01, 0.1, 1.0, 5.0, 10.0, 25.0, 50.0,
                75.0, 100.0, 125.0, 150.0,
                critical_radius * 0.999,
                critical_radius,
                critical_radius * 1.001
            };

            for (double amplitude : amplitudes) {
                for (double sign : {-1.0, 1.0}) {
                    const double h = sign * amplitude;

                    const auto nonlinear =
                        run_forward_scheduled(
                            forecast_cfg, perturb(h), future);

                    if (nonlinear.info != 0) {
                        std::cout
                            << "  e59b_h=" << h
                            << " status=NONLINEAR_SOLVE_FAILED"
                            << " info=" << nonlinear.info << '\n';
                        continue;
                    }

                    require(nonlinear.outlet_pressure.size() ==
                            nominal_forward.outlet_pressure.size(),
                            "E5.9b forecast size mismatch");

                    double max_remainder = 0.0;
                    double max_linear_change = 0.0;

                    double min_actual_pressure =
                        std::numeric_limits<double>::infinity();

                    double min_linear_margin =
                        std::numeric_limits<double>::infinity();

                    double min_nonlinear_margin =
                        std::numeric_limits<double>::infinity();

                    for (Eigen::Index k = 0;
                         k < predicted.size(); ++k) {

                        const auto j =
                            static_cast<std::size_t>(k);

                        const double nominal =
                            nominal_forward.outlet_pressure[j];

                        const double actual =
                            nonlinear.outlet_pressure[j];

                        const double linear =
                            nominal + h * predicted[k];

                        const double remainder =
                            actual - linear;

                        const double measurement_allowance =
                            measurement_radius
                            * d.measurement_unit_radius[j];

                        max_remainder =
                            std::max(max_remainder,
                                     std::abs(remainder));

                        max_linear_change =
                            std::max(max_linear_change,
                                     std::abs(h * predicted[k]));

                        min_actual_pressure =
                            std::min(min_actual_pressure, actual);

                        min_linear_margin =
                            std::min(min_linear_margin,
                                     linear
                                     - measurement_allowance
                                     - pressure_limit);

                        min_nonlinear_margin =
                            std::min(min_nonlinear_margin,
                                     actual
                                     - measurement_allowance
                                     - pressure_limit);
                    }

                    const double relative_remainder =
                        max_linear_change > 0.0
                        ? max_remainder / max_linear_change
                        : 0.0;

                    std::cout
                        << "  e59b_h=" << h
                        << " status=OK"
                        << " max_remainder_pa="
                        << max_remainder
                        << " max_relative_remainder="
                        << relative_remainder
                        << " min_actual_pressure_pa="
                        << min_actual_pressure
                        << " min_linear_margin_pa="
                        << min_linear_margin
                        << " min_nonlinear_margin_pa="
                        << min_nonlinear_margin
                        << '\n';
                }
            }



            // E5.9d: domain-aware, sampled directional crossing search.
            // This is not a nonlinear uncertainty-set certificate.
            {
                const std::vector<double> lead_times = {
                    0.0, 1.0, 2.0, 3.0, 5.0
                };

                constexpr double maximum_radius = 1000.0;
                constexpr int scan_intervals = 256;
                constexpr int bisection_iterations = 60;
                constexpr double radius_tolerance = 1.0e-6;

                enum class EvaluationStatus {
                    valid,
                    invalid_initial_state,
                    solver_failure
                };

                struct Evaluation {
                    EvaluationStatus status;
                    double margin;
                    std::size_t active_index;
                };

                auto evaluate = [&](double h, double lead) {
                    const auto initial = perturb(h);

                    for (double pressure : initial.pressure) {
                        if (!std::isfinite(pressure) ||
                            pressure <= 0.0) {
                            return Evaluation{
                                EvaluationStatus::invalid_initial_state,
                                std::numeric_limits<double>::quiet_NaN(),
                                d.time.size()
                            };
                        }
                    }

                    for (double flow : initial.flow) {
                        if (!std::isfinite(flow)) {
                            return Evaluation{
                                EvaluationStatus::invalid_initial_state,
                                std::numeric_limits<double>::quiet_NaN(),
                                d.time.size()
                            };
                        }
                    }

                    if (!std::isfinite(initial.inlet_flow)) {
                        return Evaluation{
                            EvaluationStatus::invalid_initial_state,
                            std::numeric_limits<double>::quiet_NaN(),
                            d.time.size()
                        };
                    }

                    try {
                        const auto forecast = run_forward_scheduled(
                            forecast_cfg, initial, future);

                        if (forecast.info != 0 ||
                            forecast.outlet_pressure.size() !=
                                d.time.size()) {
                            return Evaluation{
                                EvaluationStatus::solver_failure,
                                std::numeric_limits<double>::quiet_NaN(),
                                d.time.size()
                            };
                        }

                        double minimum =
                            std::numeric_limits<double>::infinity();

                        std::size_t active = d.time.size();

                        for (std::size_t k = 0;
                             k < d.time.size(); ++k) {

                            if (d.time[k] < lead)
                                continue;

                            const double pressure =
                                forecast.outlet_pressure[k];

                            if (!std::isfinite(pressure)) {
                                return Evaluation{
                                    EvaluationStatus::solver_failure,
                                    std::numeric_limits<double>::quiet_NaN(),
                                    d.time.size()
                                };
                            }

                            const double allowance =
                                measurement_radius *
                                d.measurement_unit_radius[k];

                            const double margin =
                                pressure - allowance - pressure_limit;

                            if (margin < minimum) {
                                minimum = margin;
                                active = k;
                            }
                        }

                        if (!std::isfinite(minimum)) {
                            return Evaluation{
                                EvaluationStatus::solver_failure,
                                std::numeric_limits<double>::quiet_NaN(),
                                d.time.size()
                            };
                        }

                        return Evaluation{
                            EvaluationStatus::valid,
                            minimum,
                            active
                        };

                    } catch (const std::exception&) {
                        return Evaluation{
                            EvaluationStatus::solver_failure,
                            std::numeric_limits<double>::quiet_NaN(),
                            d.time.size()
                        };
                    }
                };

                for (double lead : lead_times) {
                    for (double sign : {-1.0, 1.0}) {

                        double linear_radius =
                            std::numeric_limits<double>::infinity();

                        std::size_t linear_active = d.time.size();

                        for (std::size_t k = 0;
                             k < d.time.size(); ++k) {

                            if (d.time[k] < lead)
                                continue;

                            const double margin =
                                d.nominal_outlet_pressure[k]
                                - measurement_radius *
                                  d.measurement_unit_radius[k]
                                - pressure_limit;

                            const double slope =
                                sign * predicted[
                                    static_cast<Eigen::Index>(k)];

                            if (slope < 0.0) {
                                const double radius =
                                    margin / (-slope);

                                if (radius < linear_radius) {
                                    linear_radius = radius;
                                    linear_active = k;
                                }
                            }
                        }

                        const auto nominal = evaluate(0.0, lead);

                        require(
                            nominal.status == EvaluationStatus::valid &&
                            nominal.margin > 0.0,
                            "E5.9d invalid nominal decision");

                        double previous_radius = 0.0;
                        double previous_margin = nominal.margin;

                        double crossing_lower = 0.0;
                        double crossing_upper = 0.0;

                        bool crossed = false;
                        bool domain_limit = false;
                        bool solver_failure = false;

                        double last_valid_radius = 0.0;
                        std::size_t crossing_active = d.time.size();

                        for (int i = 1; i <= scan_intervals; ++i) {
                            const double radius =
                                maximum_radius *
                                static_cast<double>(i) /
                                static_cast<double>(scan_intervals);

                            const auto result =
                                evaluate(sign * radius, lead);

                            if (result.status ==
                                EvaluationStatus::invalid_initial_state) {
                                domain_limit = true;
                                break;
                            }

                            if (result.status ==
                                EvaluationStatus::solver_failure) {
                                solver_failure = true;
                                break;
                            }

                            last_valid_radius = radius;

                            if (previous_margin > 0.0 &&
                                result.margin <= 0.0) {
                                crossing_lower = previous_radius;
                                crossing_upper = radius;
                                crossing_active = result.active_index;
                                crossed = true;
                                break;
                            }

                            previous_radius = radius;
                            previous_margin = result.margin;
                        }

                        bool refinement_failed = false;

                        if (crossed) {
                            for (int iteration = 0;
                                 iteration < bisection_iterations;
                                 ++iteration) {

                                const double middle =
                                    0.5 * (crossing_lower +
                                           crossing_upper);

                                const auto result =
                                    evaluate(sign * middle, lead);

                                if (result.status !=
                                    EvaluationStatus::valid) {
                                    refinement_failed = true;
                                    break;
                                }

                                if (result.margin > 0.0) {
                                    crossing_lower = middle;
                                } else {
                                    crossing_upper = middle;
                                    crossing_active =
                                        result.active_index;
                                }

                                if (crossing_upper -
                                    crossing_lower <
                                    radius_tolerance) {
                                    break;
                                }
                            }
                        }

                        const bool valid_crossing =
                            crossed && !refinement_failed;

                        const double nonlinear_radius =
                            valid_crossing
                            ? 0.5 * (crossing_lower +
                                     crossing_upper)
                            : std::numeric_limits<double>::quiet_NaN();

                        const char* status = "RANGE_EXHAUSTED";

                        if (refinement_failed)
                            status = "REFINEMENT_FAILURE";
                        else if (valid_crossing)
                            status = "CROSSING";
                        else if (domain_limit)
                            status = "DOMAIN_LIMIT";
                        else if (solver_failure)
                            status = "SOLVER_FAILURE";

                        std::cout
                            << "  e59d_lead=" << lead
                            << " sign=" << sign
                            << " status=" << status
                            << " linear_radius=" << linear_radius
                            << " linear_active_index="
                            << linear_active
                            << " nonlinear_radius="
                            << nonlinear_radius
                            << " nonlinear_active_index="
                            << crossing_active
                            << " last_valid_radius="
                            << last_valid_radius
                            << " radius_difference="
                            << (valid_crossing &&
                                std::isfinite(linear_radius)
                                ? nonlinear_radius - linear_radius
                                : std::numeric_limits<double>::quiet_NaN())
                            << '\n';
                    }
                }

                std::cout
                    << "  e59d_lead_time_analysis=COMPLETE\n";
            }

            std::cout
                << "  e59b_finite_amplitude_sweep=COMPLETE\n";
        }

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
