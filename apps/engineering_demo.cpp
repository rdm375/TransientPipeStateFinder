#include "tpsf/compiled_provenance.hpp"

#include "tpsf/dense_full_space_sqp.hpp"
#include "tpsf/decision_timescales.hpp"
#include "tpsf/engineering_report.hpp"
#include "tpsf/state_qualification.hpp"
#include "tpsf/decision_boundary_sweep.hpp"
#include "tpsf/forecast_attenuation.hpp"
#include "tpsf/sensitivity_decay.hpp"
#include "tpsf/physical_boundary_channels.hpp"
#include "tpsf/forecast_continuation.hpp"
#include "tpsf/pipeline_decision.hpp"
#include "tpsf/scheduled_forward.hpp"
#include "tpsf/statefinding_quality.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
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

Configuration configuration()
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
    require(
        forward.summary.info == 0 &&
        forward.summary.completed_steps == cfg.steps,
        "Forward integration failed");

    FullSpaceTrajectory trajectory;
    trajectory.n = cfg.cells;

    for (const auto& snapshot :
         forward.trajectory.snapshots) {

        trajectory.states.push_back(snapshot.state);

        trajectory.inputs.push_back(
            snapshot.inputs);

        trajectory.reconstructed.push_back(
            snapshot.physical_channels);
    }

    require(
        trajectory.states.size() ==
            static_cast<std::size_t>(cfg.steps + 1),
        "Incomplete forward trajectory");

    return trajectory;
}

FullSpaceProblem make_problem(
    const FullSpaceTrajectory& reference,
    const FullSpaceTrajectory& truth,
    const pipe_sim::TransientParameters& parameters)
{
    FullSpaceProblemConfiguration cfg;

    cfg.reference = reference;
    cfg.parameters = parameters;

    const std::size_t levels =
        reference.states.size();

    const std::size_t state_size =
        static_cast<std::size_t>(
            2 * reference.n + 2);

    const std::size_t block_size =
        state_size + 6;

    cfg.variables.assign(
        levels * block_size,
        ProblemVariable{});

    for (std::size_t k = 0; k < levels; ++k) {

        const std::size_t offset =
            k * block_size;

        for (std::size_t j = 0;
             j < state_size; ++j) {

            auto& variable =
                cfg.variables[offset + j];

            // Two independently unknown initial
            // pressures; all later states are free.
            variable.free =
                k != 0 || j == 3 || j == 5;

            const bool pressure =
                j == 0 ||
                (j >= 3 && j % 2 == 1);

            variable.scale =
                pressure ? 1.0e6 : 10.0;
        }

        // Prescribed controls remain fixed.
        cfg.variables[
            offset + state_size].free = false;

        cfg.variables[
            offset + state_size + 1].free = false;

        // Four independently reconstructed
        // physical boundary channels.
        for (std::size_t channel = 0;
             channel < 4; ++channel) {

            auto& variable =
                cfg.variables[
                    offset +
                    state_size +
                    2 +
                    channel];

            variable.free = true;

            variable.scale =
                channel == 0 || channel == 2
                    ? 1.0e6
                    : 10.0;
        }

        BoundaryMeasurement measurement;

        measurement.observed =
            truth.reconstructed[k];

        // R1: All four observation channels.
        measurement.available = {
            true, true, true, true
        };

        measurement.sigma = {
            1.0e3,
            0.1,
            1.0e3,
            0.1
        };

        cfg.measurements.push_back(
            measurement);
    }

    const auto constraints =
        assemble_full_space_trajectory_constraints(
            reference,
            parameters);

    cfg.constraint_scales.assign(
        constraints.combined().size(),
        1.0);

    return FullSpaceProblem(
        std::move(cfg));
}

InitialState perturb_initial_state(
    InitialState state,
    std::size_t coordinate,
    double displacement)
{
    switch (coordinate) {

    case 0:
        state.inlet_flow += displacement;
        break;

    case 1:
        state.flow[0] += displacement;
        break;

    case 2:
        state.pressure[1] += displacement;
        break;

    case 3:
        state.flow[1] += displacement;
        break;

    case 4:
        state.pressure[2] += displacement;
        break;

    default:
        throw std::runtime_error(
            "Invalid initial-state coordinate");
    }

    return state;
}

struct MemorySample {
    double time = 0.0;
    double state_gain = 0.0;
    double boundary_gain = 0.0;
    double outlet_pressure_gain = 0.0;
};

struct MemoryResult {
    std::vector<MemorySample> samples;
    double maximum_tangent_residual = 0.0;
    double maximum_tangent_mass_defect = 0.0;
};


MemoryResult analyze_initial_state_memory(
    const Configuration& cfg,
    const InitialState& nominal,
    const BoundarySchedule& future)
{
    const auto baseline =
        recorded_trajectory(
            cfg,
            run_forward_scheduled_recorded(
                cfg, nominal, future));

    constexpr std::array<int, 5> indices{
        1, 2, 3, 4, 5
    };

    constexpr std::array<double, 5> scales{
        10.0, 10.0, 1.0e6, 10.0, 1.0e6
    };

    MemoryResult result;
    result.samples.resize(baseline.states.size());

    std::array<ForecastTangentHistory, 5> histories;

    for (std::size_t j = 0; j < indices.size(); ++j) {

        Eigen::VectorXd direction =
            Eigen::VectorXd::Zero(6);

        direction[indices[j]] = scales[j];

        histories[j] =
            propagate_forecast_tangent(
                baseline,
                direction,
                cfg.parameters);

        require(
            histories[j].states.size() ==
                result.samples.size(),
            "Incomplete tangent history");

        result.maximum_tangent_residual =
            std::max(
                result.maximum_tangent_residual,
                histories[j].maximum_equation_residual);

        result.maximum_tangent_mass_defect =
            std::max(
                result.maximum_tangent_mass_defect,
                histories[j].maximum_mass_balance_defect);
    }

    for (std::size_t k = 0;
         k < result.samples.size(); ++k) {

        Eigen::MatrixXd S(6, 5);
        Eigen::MatrixXd H(4, 5);

        for (std::size_t j = 0;
             j < indices.size(); ++j) {

            const auto& tangent =
                histories[j].states[k];

            std::array<double, 6> state_tangent{};

            for (int i = 0; i < 6; ++i) {

                const bool pressure =
                    i == 0 ||
                    (i >= 3 && i % 2 == 1);

                const double output_scale =
                    pressure ? 1.0e6 : 10.0;

                S(i, j) =
                    tangent[i] / output_scale;

                state_tangent[i] = tangent[i];
            }

            const PhysicalBoundaryInputs zero_inputs{};

            const auto channels =
                physical_boundary_channels_jvp(
                    cfg.cells,
                    state_tangent,
                    zero_inputs);

            for (int channel = 0;
                 channel < 4; ++channel) {

                const double output_scale =
                    channel == 0 || channel == 2
                        ? 1.0e6 : 10.0;

                H(channel, j) =
                    channels[channel] /
                    output_scale;
            }
        }

        auto& sample = result.samples[k];

        sample.time = k * cfg.parameters.dt;

        sample.state_gain =
            Eigen::JacobiSVD<Eigen::MatrixXd>(S)
                .singularValues()[0];

        sample.boundary_gain =
            Eigen::JacobiSVD<Eigen::MatrixXd>(H)
                .singularValues()[0];

        sample.outlet_pressure_gain =
            H.row(2).norm();
    }

    return result;
}


void write_reports(
    const std::filesystem::path& directory,
    const DenseSQPResult& solution,
    const StatefindingQuality& quality,
    const PipelineDecisionResult& decision,
    const MemoryResult& memory)
{

    std::vector<double> times;
    std::vector<double> state_gains;
    std::vector<double> boundary_gains;
    std::vector<double> outlet_gains;

    for (const auto& sample : memory.samples) {

        times.push_back(sample.time);
        state_gains.push_back(sample.state_gain);
        boundary_gains.push_back(sample.boundary_gain);
        outlet_gains.push_back(
            sample.outlet_pressure_gain);
    }

    SensitivityDecayOptions options;
    options.relative_threshold = 0.05;

    const auto state_decay =
        assess_sensitivity_decay(
            times, state_gains, options);

    const auto boundary_decay =
        assess_sensitivity_decay(
            times, boundary_gains, options);

    const auto outlet_decay =
        assess_sensitivity_decay(
            times, outlet_gains, options);

    std::filesystem::create_directories(
        directory);

    std::ofstream json(
        directory / "engineering_demo.json");

    std::ofstream markdown(
        directory / "engineering_demo.md");

    require(
        json.good() && markdown.good(),
        "Unable to create engineering reports");

    json << std::setprecision(17);

    json
        << "{\n"
        << "  \"schema\": "
           "\"tpsf.engineering-demo.v0.1\",\n"
        << "  \"synthetic_data\": true,\n"
        << "  \"forward_solver_revision\": "
           "\"670ac493e4461e678ca66f2254dc3d402df3f312\",\n"
        << "  \"four_channel_observations\": true,\n"
        << "  \"sqp_converged\": "
        << (solution.converged()
                ? "true" : "false")
        << ",\n"
        << "  \"constraint_l2\": "
        << solution.final_constraint_l2
        << ",\n"
        << "  \"feasible_dimension\": "
        << quality.feasible_dimension
        << ",\n"
        << "  \"measurement_rank\": "
        << quality.measurement_rank
        << ",\n"
        << "  \"unresolved_directions\": "
        << quality.unresolved_directions
        << ",\n"
        << "  \"decision_classification\": "
        << static_cast<int>(
            decision.decision.classification)
        << ",\n"
        << "  \"memory_tangent_residual\": "
        << memory.maximum_tangent_residual
        << ",\n"
        << "  \"memory_tangent_mass_defect\": "
        << memory.maximum_tangent_mass_defect
        << ",\n"
        << "  \"attenuation_relative_threshold\": 0.05,\n"
        << "  \"state_persistent_attenuation_time_s\": "
        << (state_decay.relative.time
                ? std::to_string(*state_decay.relative.time)
                : "null")
        << ",\n"
        << "  \"boundary_persistent_attenuation_time_s\": "
        << (boundary_decay.relative.time
                ? std::to_string(*boundary_decay.relative.time)
                : "null")
        << ",\n"
        << "  \"outlet_persistent_attenuation_time_s\": "
        << (outlet_decay.relative.time
                ? std::to_string(*outlet_decay.relative.time)
                : "null")
        << ",\n"
        << "  \"memory\": [\n";

    for (std::size_t k = 0;
         k < memory.samples.size(); ++k) {

        const auto& sample =
            memory.samples[k];

        json
            << "    {\n"
            << "      \"time_s\": "
            << sample.time << ",\n"
            << "      \"state_gain\": "
            << sample.state_gain << ",\n"
            << "      \"boundary_gain\": "
            << sample.boundary_gain << ",\n"
            << "      \"outlet_pressure_gain\": "
            << sample.outlet_pressure_gain
            << "\n"
            << "    }";

        if (k + 1 < memory.samples.size())
            json << ",";

        json << "\n";
    }

    json << "  ],\n";

    const auto& directional =
        decision.directional;

    json << "  \"forecast\": [\n";

    for (std::size_t k = 0;
         k < directional.time.size(); ++k) {

        json
            << "    {\n"
            << "      \"time_s\": "
            << directional.time[k] << ",\n"
            << "      \"outlet_pressure_pa\": "
            << directional.nominal_outlet_pressure[k]
            << ",\n"
            << "      \"measurement_unit_radius_pa\": "
            << directional.measurement_unit_radius[k]
            << ",\n"
            << "      \"unresolved_unit_radius_pa\": "
            << directional.unresolved_unit_radius[k]
            << "\n"
            << "    }";

        if (k + 1 < directional.time.size())
            json << ",";

        json << "\n";
    }

    json << "  ]\n}\n";

    markdown
        << "# TransientPipeStateFinder\n\n"
        << "## Synthetic engineering demonstration\n\n"
        << "**Status:** Research demonstration. "
           "Not an operational certificate.\n\n"
        << "### Reconstruction\n\n"
        << "- SQP converged: "
        << (solution.converged() ? "yes" : "no")
        << "\n"
        << "- Constraint L2: "
        << solution.final_constraint_l2
        << "\n"
        << "- Feasible dimension: "
        << quality.feasible_dimension
        << "\n"
        << "- Measurement rank: "
        << quality.measurement_rank
        << "\n"
        << "- Unresolved directions: "
        << quality.unresolved_directions
        << "\n\n"
        << "### Initial-state memory\n\n"
        << "Finite-difference approximation "
           "of the initial-state sensitivity operator. "
           "Pressure scale: 1 MPa. "
           "Flow scale: 10 kg/s. "
           "Future controls remain fixed.\n\n"
        << "Maximum exact tangent equation residual: "
        << memory.maximum_tangent_residual
        << "\n\n"

        << "### Persistent 5% attenuation\n\n"
        << "- Full state: "
        << (state_decay.relative.time
                ? std::to_string(*state_decay.relative.time) + " s"
                : "not attained")
        << "\n"
        << "- Four-channel boundary response: "
        << (boundary_decay.relative.time
                ? std::to_string(*boundary_decay.relative.time) + " s"
                : "not attained")
        << "\n"
        << "- Outlet pressure: "
        << (outlet_decay.relative.time
                ? std::to_string(*outlet_decay.relative.time) + " s"
                : "not attained")
        << "\n\n"
        << "| Time (s) | State gain | "
           "Boundary gain | Outlet pressure gain |\n"
        << "|---:|---:|---:|---:|\n";

    for (const auto& sample :
         memory.samples) {

        markdown
            << "| " << sample.time
            << " | " << sample.state_gain
            << " | " << sample.boundary_gain
            << " | " << sample.outlet_pressure_gain
            << " |\n";
    }

    markdown
        << "\n### Engineering decision\n\n"
        << "Minimum outlet pressure: 4 MPa.\n\n"
        << "Classification: "
        << static_cast<int>(
            decision.decision.classification)
        << " (0=satisfactory, "
           "1=violation, 2=indeterminate).\n\n"
        << "Measurement radius is a deterministic "
           "bound on standardized perturbations, "
           "not a statistical confidence level.\n\n"
        << "Unresolved directions are unbounded.\n\n"
        << "Classification applies to the "
           "local-linearized sampled-time model. "
           "It is not a nonlinear certificate.\n\n"
        << "### Limitations\n\n"
        << "Synthetic two-cell example. "
           "The inlet pressure and outlet flow "
           "are prescribed controls and do not "
           "independently identify the initial state. "
           "The memory operator is estimated "
           "using nonlinear finite differences, "
           "not the exact discrete tangent-linear "
           "propagator. Persistent attenuation times "
           "and grid convergence remain future work.\n";
}


// R5.1: causal state reconstruction and independent
// nonlinear restart qualification.
//
// At analysis step k, only measurements at levels 0..k
// are supplied to the statefinder.
//
// The synthetic truth is used solely for validation,
// never for the initial guess or the statefinder's
// measurement history beyond the candidate time.

StateQualificationResult run_qualification_experiment()
{
    auto cfg = configuration();

    constexpr int analysis_steps = 8;
    constexpr int restart_steps = 4;
    constexpr int total_steps =
        analysis_steps + restart_steps;

    cfg.steps = total_steps;

    const InitialState truth_initial{
        {5.0e6, 5.0e6, 5.0e6},
        {10.0, 10.0},
        10.0
    };

    InitialState guessed_initial = truth_initial;

    guessed_initial.pressure[1] += 1.5e4;
    guessed_initial.pressure[2] += 2.0e4;

    BoundarySchedule schedule;

    for (int k = 0; k <= total_steps; ++k) {
        schedule.inlet_pressure.push_back(5.0e6);
        schedule.outlet_flow.push_back(
            10.0 + 0.02 * k);
    }

    const auto truth =
        recorded_trajectory(
            cfg,
            run_forward_scheduled_recorded(
                cfg,
                truth_initial,
                schedule));

    const auto guess =
        recorded_trajectory(
            cfg,
            run_forward_scheduled_recorded(
                cfg,
                guessed_initial,
                schedule));

    StateQualificationOptions qualification_options;

    qualification_options.pressure_tolerance_pa =
        1.0e4;

    qualification_options.flow_tolerance_kg_s =
        0.1;

    qualification_options.restart_horizon_s =
        restart_steps * cfg.parameters.dt;

    std::vector<StateQualificationSample> samples;

    for (int k = 1; k <= analysis_steps; ++k) {

        StateQualificationSample sample;

        sample.analysis_time_s =
            k * cfg.parameters.dt;

        try {
            // Only the causal measurement prefix is
            // visible to this reconstruction.
            FullSpaceTrajectory truth_prefix = truth;
            FullSpaceTrajectory guess_prefix = guess;

            const std::size_t levels =
                static_cast<std::size_t>(k + 1);

            truth_prefix.states.resize(levels);
            truth_prefix.inputs.resize(levels);
            truth_prefix.reconstructed.resize(levels);

            guess_prefix.states.resize(levels);
            guess_prefix.inputs.resize(levels);
            guess_prefix.reconstructed.resize(levels);

            const auto problem =
                make_problem(
                    guess_prefix,
                    truth_prefix,
                    cfg.parameters);

            const Eigen::VectorXd initial_coordinates =
                Eigen::VectorXd::Zero(
                    static_cast<Eigen::Index>(
                        problem.free_variable_count()));

            DenseSQPOptions options;
            options.maximum_iterations = 40;
            options.feasibility_tolerance = 1.0e-8;
            options.stationarity_tolerance = 1.0e-8;

            const auto solution =
                solve_dense_full_space_sqp(
                    problem,
                    initial_coordinates,
                    options);

            sample.reconstruction_converged =
                solution.converged();

            if (!solution.converged()) {
                samples.push_back(sample);
                continue;
            }

            const std::vector<double> coordinates(
                solution.coordinates.data(),
                solution.coordinates.data() +
                    solution.coordinates.size());

            const auto reconstructed =
                problem.trajectory(coordinates);

            const auto restart_initial =
                initial_state_from_full_space(
                    reconstructed,
                    static_cast<std::size_t>(k));

            auto restart_cfg = cfg;
            restart_cfg.steps = restart_steps;

            BoundarySchedule restart_schedule;

            for (int j = 0; j <= restart_steps; ++j) {
                restart_schedule.inlet_pressure.push_back(
                    schedule.inlet_pressure[k + j]);

                restart_schedule.outlet_flow.push_back(
                    schedule.outlet_flow[k + j]);
            }

            // Independent nonlinear integration from
            // the statefinder's reconstructed state.
            const auto restarted =
                recorded_trajectory(
                    restart_cfg,
                    run_forward_scheduled_recorded(
                        restart_cfg,
                        restart_initial,
                        restart_schedule));

            double maximum_pressure_error = 0.0;
            double maximum_flow_error = 0.0;

            for (int j = 0; j <= restart_steps; ++j) {

                const auto& estimated =
                    restarted.states[
                        static_cast<std::size_t>(j)];

                const auto& reference =
                    truth.states[
                        static_cast<std::size_t>(k + j)];

                require(
                    estimated.size() == reference.size(),
                    "Restart state size mismatch");

                // State layout for n=2:
                // [inlet pressure, inlet flow,
                //  cell flow, pressure, cell flow,
                //  outlet pressure].
                //
                // More generally, pressure indices are
                // 0,3,5,...,2*n+1.
                // Remaining entries are mass flows.
                for (std::size_t i = 0;
                     i < estimated.size(); ++i) {

                    const bool pressure =
                        i == 0 ||
                        (i >= 3 && i % 2 == 1);

                    const double error =
                        std::abs(
                            estimated[i] - reference[i]);

                    if (pressure) {
                        maximum_pressure_error =
                            std::max(
                                maximum_pressure_error,
                                error);
                    } else {
                        maximum_flow_error =
                            std::max(
                                maximum_flow_error,
                                error);
                    }
                }
            }

            sample.pressure_error_pa =
                maximum_pressure_error;

            sample.flow_error_kg_s =
                maximum_flow_error;

            sample.restart_succeeded = true;

        } catch (const std::exception& error) {

            std::cerr
                << "Qualification candidate at "
                << sample.analysis_time_s
                << " s failed: "
                << error.what() << '\n';
        }

        samples.push_back(sample);
    }

    return assess_state_qualification(
        std::move(samples),
        qualification_options);
}

} // namespace

int main(int argc, char** argv)
{
    try {

        const std::filesystem::path output =
            argc > 1
                ? argv[1]
                : "reports/engineering-demo";

        // R5.1 primary engineering result.
        const auto qualification =
            run_qualification_experiment();

        write_state_qualification(
            qualification,
            output);

        std::cout
            << "State qualification: "
            << (qualification.qualified
                ? "QUALIFIED"
                : "NOT QUALIFIED")
            << '\n';

        if (qualification.qualification_time_s) {
            std::cout
                << "Qualification time: "
                << *qualification.qualification_time_s
                << " s\n";
        }

        const auto cfg =
            configuration();

        const InitialState initial{
            {5.0e6, 5.0e6, 5.0e6},
            {10.0, 10.0},
            10.0
        };

        const BoundarySchedule schedule{
            {5.0e6, 5.0e6, 5.0e6},
            {10.0, 10.1, 10.2}
        };

        const auto truth =
            recorded_trajectory(
                cfg,
                run_forward_scheduled_recorded(
                    cfg, initial, schedule));

        InitialState guessed_initial =
            initial;

        guessed_initial.pressure[1] +=
            1.5e4;

        guessed_initial.pressure[2] +=
            2.0e4;

        const auto reference =
            recorded_trajectory(
                cfg,
                run_forward_scheduled_recorded(
                    cfg,
                    guessed_initial,
                    schedule));

        const auto problem =
            make_problem(
                reference,
                truth,
                cfg.parameters);

        const Eigen::VectorXd initial_coordinates =
            Eigen::VectorXd::Zero(
                static_cast<Eigen::Index>(
                    problem.free_variable_count()));

        DenseSQPOptions options;

        options.maximum_iterations = 40;
        options.feasibility_tolerance = 1.0e-8;
        options.stationarity_tolerance = 1.0e-8;

        const auto solution =
            solve_dense_full_space_sqp(
                problem,
                initial_coordinates,
                options);

        require(
            solution.converged(),
            "SQP reconstruction did not converge");

        const auto quality =
            assess_statefinding_quality(
                problem,
                solution);

        const std::vector<double> coordinates(
            solution.coordinates.data(),
            solution.coordinates.data() +
                solution.coordinates.size());

        const auto reconstructed =
            problem.trajectory(
                coordinates);

        const auto forecast_initial =
            initial_state_from_full_space(
                reconstructed,
                reconstructed.states.size() - 1);

        auto forecast_configuration =
            cfg;

        forecast_configuration.steps = 10;

        BoundarySchedule future;

        for (int k = 0;
             k <= forecast_configuration.steps;
             ++k) {

            const double time =
                k * cfg.parameters.dt;

            future.inlet_pressure.push_back(
                5.0e6);

            future.outlet_flow.push_back(
                10.2 + 0.02 * time);
        }

        PipelineDecisionOptions decision_options;

        decision_options.minimum_outlet_pressure_pa =
            4.0e6;

        decision_options.measurement_radius =
            1.0;

        decision_options.bound_unresolved =
            false;

        const auto decision =
            assess_pipeline_decision(
                problem,
                solution.coordinates,
                quality,
                forecast_configuration,
                future,
                decision_options);

        const auto timescales =
            assess_operational_timescales(
                decision);

        const auto memory =
            analyze_initial_state_memory(
                forecast_configuration,
                forecast_initial,
                future);

        write_reports(
            output,
            solution,
            quality,
            decision,
            memory);

        // R3.1 reusable engineering reporting.
        EngineeringReport assessment;

        assessment.quality = quality;
        assessment.decision = decision;

        assessment.scenario.cells = cfg.cells;
        assessment.scenario.reconstruction_steps = cfg.steps;
        assessment.scenario.forecast_steps =
            forecast_configuration.steps;

        assessment.scenario.dt = cfg.parameters.dt;
        assessment.scenario.theta = cfg.parameters.theta;

        assessment.scenario.minimum_outlet_pressure_pa =
            decision_options.minimum_outlet_pressure_pa;

        assessment.scenario.measurement_radius =
            decision_options.measurement_radius;

        assessment.scenario.bound_unresolved =
            decision_options.bound_unresolved;

        assessment.scenario.unresolved_radius =
            decision_options.unresolved_radius;

        // The synthetic demonstration has four available
        // channels at every reconstruction time level.
        for (std::size_t k = 0;
             k < truth.states.size(); ++k) {

            EngineeringMeasurement measurement;

            measurement.available = {
                true, true, true, true
            };

            measurement.sigma = {
                1.0e3, 0.1, 1.0e3, 0.1
            };

            assessment.scenario.measurements.push_back(
                measurement);
        }

        assessment.maximum_memory_tangent_residual =
            memory.maximum_tangent_residual;

        assessment.maximum_memory_mass_defect =
            memory.maximum_tangent_mass_defect;

        std::vector<double> memory_times;
        std::vector<double> state_gains;
        std::vector<double> boundary_gains;
        std::vector<double> outlet_gains;

        for (const auto& sample : memory.samples) {
            memory_times.push_back(sample.time);
            state_gains.push_back(sample.state_gain);
            boundary_gains.push_back(sample.boundary_gain);
            outlet_gains.push_back(
                sample.outlet_pressure_gain);
        }

        SensitivityDecayOptions attenuation_options;
        attenuation_options.relative_threshold = 0.05;

        const auto state_attenuation =
            assess_sensitivity_decay(
                memory_times,
                state_gains,
                attenuation_options);

        const auto boundary_attenuation =
            assess_sensitivity_decay(
                memory_times,
                boundary_gains,
                attenuation_options);

        const auto outlet_attenuation =
            assess_sensitivity_decay(
                memory_times,
                outlet_gains,
                attenuation_options);

        assessment.state_attenuation_time =
            state_attenuation.relative.time;

        assessment.boundary_attenuation_time =
            boundary_attenuation.relative.time;

        assessment.outlet_attenuation_time =
            outlet_attenuation.relative.time;

        // Provenance was embedded by the build, not queried at report time.
        assessment.provenance = compiled_engineering_provenance();

        assessment.qualification = qualification;
        write_engineering_report(
            assessment,
            output);

        // R4.1: reuse one nonlinear forecast and its
        // local-linearized uncertainty bounds.

        const std::vector<double> thresholds{
            4.0e6,
            4.9e6,
            4.99e6,
            4.999e6,
            5.0e6,
            5.001e6,
            5.01e6
        };

        const auto boundary_sweep =
            assess_decision_boundary_sweep(
                decision,
                thresholds);

        write_decision_boundary_sweep(
            boundary_sweep,
            output);




        std::cout
            << "Engineering reports written to "
            << output << '\n'
            << "Measurement rank: "
            << quality.measurement_rank << '\n'
            << "Unresolved directions: "
            << quality.unresolved_directions << '\n'
            << "Memory tangent residual: "
            << memory.maximum_tangent_residual
            << '\n';

        (void)timescales;

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << "Engineering demonstration failed: "
            << error.what() << '\n';

        return 1;
    }
}
