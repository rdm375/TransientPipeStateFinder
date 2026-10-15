#include "tpsf/full_space_problem.hpp"
#include "tpsf/physical_boundary_channels.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
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

void near(double a, double b, double atol,
          double rtol, const std::string& message)
{
    require(std::isfinite(a) && std::isfinite(b),
            message + ": nonfinite");

    if (std::abs(a-b) >
        atol + rtol*std::max(std::abs(a), std::abs(b))) {
        throw std::runtime_error(
            message + ": " + std::to_string(a) +
            " vs " + std::to_string(b));
    }
}

template<class F>
void expect_invalid(F&& f)
{
    bool caught = false;

    try {
        f();
    }
    catch (const std::invalid_argument&) {
        caught = true;
    }

    require(caught, "Expected invalid_argument");
}

double dot(const std::vector<double>& a,
           const std::vector<double>& b)
{
    require(a.size() == b.size(), "Dot size mismatch");

    double result = 0.0;

    for (std::size_t i = 0; i < a.size(); ++i)
        result += a[i]*b[i];

    return result;
}

FullSpaceProblemConfiguration configuration(
    int n = 5,
    int steps = 3)
{
    FullSpaceProblemConfiguration cfg;

    cfg.parameters = pipe_sim::TransientParameters{
        100.0, 1.0, 0.65,
        0.5, 0.19634954084936207,
        288.15, 0.9, 518.3,
        1.1e-5, 1.0e-5
    };

    cfg.reference.n = n;

    const int nu = 2*n+2;
    const int block = nu+6;

    for (int k = 0; k <= steps; ++k) {
        std::vector<double> state(nu);

        for (int i = 0; i < nu; ++i) {
            const bool pressure =
                i == 0 || (i >= 3 && i%2 == 1);

            state[i] = pressure
                ? 5.0e6 - 1500.0*i + 300.0*k
                : 10.0 + 0.04*i + 0.02*k;
        }

        cfg.reference.states.push_back(state);

        const PhysicalBoundaryInputs q{
            5.0e6 + 500.0*k,
            10.0 + 0.1*k
        };

        cfg.reference.inputs.push_back(q);

        auto channels =
            extract_physical_boundary_channels(
                n, state, q);

        for (std::size_t c = 0; c < 4; ++c)
            channels[c] += (c%2 == 0 ? 100.0 : 0.1);

        cfg.reference.reconstructed.push_back(channels);

        BoundaryMeasurement measurement;

        for (std::size_t c = 0; c < 4; ++c) {
            const double sigma =
                c%2 == 0 ? 200.0 : 0.2;

            measurement.observed[c] =
                channels[c] - 0.5*sigma;

            measurement.sigma[c] = sigma;

            measurement.available[c] =
                (k+c)%3 != 0;
        }

        cfg.measurements.push_back(measurement);
    }

    const std::size_t variable_count =
        static_cast<std::size_t>(steps+1)*block;

    cfg.variables.resize(variable_count);

    for (int k = 0; k <= steps; ++k) {
        const std::size_t offset =
            static_cast<std::size_t>(k)*block;

        for (int i = 0; i < nu; ++i) {
            const bool pressure =
                i == 0 || (i >= 3 && i%2 == 1);

            cfg.variables[offset+i].scale =
                pressure ? 1.0e5 : 5.0;
        }

        // Both boundary inputs are independently
        // configurable at each time level.
        cfg.variables[offset+nu].scale = 1.0e5;
        cfg.variables[offset+nu+1].scale = 5.0;

        for (int c = 0; c < 4; ++c)
            cfg.variables[offset+nu+2+c].scale =
                c%2 == 0 ? 1.0e5 : 5.0;

        // Fix the inlet pressure at every time level.
        cfg.variables[offset+nu].free = false;

        // Fix outlet flow at even time levels only.
        cfg.variables[offset+nu+1].free = k%2 != 0;
    }

    const std::size_t constraint_count =
        static_cast<std::size_t>(steps)*nu +
        4*static_cast<std::size_t>(steps+1);

    cfg.constraint_scales.resize(constraint_count);

    for (std::size_t i = 0;
         i < constraint_count; ++i) {
        cfg.constraint_scales[i] =
            i%3 == 0 ? 1.0e-5 :
            i%3 == 1 ? 0.1 : 2.0;
    }

    return cfg;
}

std::vector<double> direction(std::size_t n)
{
    std::vector<double> result(n);

    for (std::size_t i = 0; i < n; ++i)
        result[i] = std::sin(0.31*(i+1));

    return result;
}

std::vector<double> dual(std::size_t n)
{
    std::vector<double> result(n);

    for (std::size_t i = 0; i < n; ++i)
        result[i] = std::cos(0.23*(i+1));

    return result;
}

void test_round_trip()
{
    const FullSpaceProblem problem(configuration());

    const auto u = direction(problem.free_variable_count());

    const auto trajectory = problem.trajectory(u);
    const auto recovered = problem.coordinates(trajectory);

    require(recovered.size() == u.size(),
            "Round-trip dimension mismatch");

    for (std::size_t i = 0; i < u.size(); ++i)
        near(recovered[i], u[i],
             2.0e-10, 2.0e-10,
             "Scaled coordinate round trip");

    const auto reference =
        problem.trajectory(
            std::vector<double>(
                problem.free_variable_count(), 0.0));

    for (std::size_t k = 0;
         k < reference.inputs.size(); ++k) {
        require(
            trajectory.inputs[k].inlet_pressure ==
            reference.inputs[k].inlet_pressure,
            "Fixed inlet pressure changed");

        if (k%2 == 0)
            require(
                trajectory.inputs[k].outlet_flow ==
                reference.inputs[k].outlet_flow,
                "Fixed outlet flow changed");
    }

    auto changed = trajectory;
    changed.inputs[0].inlet_pressure += 1.0;

    expect_invalid([&] {
        problem.coordinates(changed);
    });
}

void test_reference_evaluation()
{
    const auto cfg = configuration();
    const FullSpaceProblem problem(cfg);

    const std::vector<double> zero(
        problem.free_variable_count(), 0.0);

    const auto result = problem.evaluate(zero);

    const auto direct_measurement =
        evaluate_full_space_measurement_objective(
            cfg.reference.reconstructed,
            cfg.measurements);

    near(result.measurement.value,
         direct_measurement.value,
         0.0, 0.0,
         "Reference objective mismatch");

    const auto direct_constraints =
        assemble_full_space_trajectory_constraints(
            cfg.reference, cfg.parameters).combined();

    require(result.constraints.size() ==
            direct_constraints.size(),
            "Constraint dimension mismatch");

    for (std::size_t i = 0;
         i < direct_constraints.size(); ++i)
        near(result.constraints[i],
             direct_constraints[i]*
                 cfg.constraint_scales[i],
             0.0, 0.0,
             "Scaled constraint mismatch");
}

void test_transpose(int n, int steps)
{
    const FullSpaceProblem problem(
        configuration(n, steps));

    auto u = direction(problem.free_variable_count());

    for (double& value : u)
        value *= 0.01;

    const auto v = direction(problem.free_variable_count());
    const auto w = dual(problem.constraint_count());

    const auto jv = problem.constraint_jvp(u, v);
    const auto jtw = problem.constraint_vjp(u, w);

    const double lhs = dot(jv, w);
    const double rhs = dot(v, jtw);

    near(lhs, rhs, 1.0e-8, 2.0e-12,
         "Scaled JVP/VJP transpose identity");
}

void test_objective_gradient()
{
    const FullSpaceProblem problem(configuration());

    auto u = direction(problem.free_variable_count());

    for (double& value : u)
        value *= 0.01;

    const auto v = direction(problem.free_variable_count());

    const auto evaluation = problem.evaluate(u);

    const double analytic =
        dot(evaluation.objective_gradient, v);

    constexpr double h = 1.0e-5;

    auto plus = u;
    auto minus = u;

    for (std::size_t i = 0; i < u.size(); ++i) {
        plus[i] += h*v[i];
        minus[i] -= h*v[i];
    }

    const double numerical =
        (problem.evaluate(plus).measurement.value -
         problem.evaluate(minus).measurement.value) /
        (2.0*h);

    near(numerical, analytic,
         2.0e-5, 2.0e-5,
         "Scaled objective gradient");
}

void test_constraint_finite_difference()
{
    const FullSpaceProblem problem(configuration());

    auto u = direction(problem.free_variable_count());

    for (double& value : u)
        value *= 0.01;

    const auto v = direction(problem.free_variable_count());

    const auto analytic = problem.constraint_jvp(u, v);

    constexpr double h = 1.0e-4;

    auto plus = u;
    auto minus = u;

    for (std::size_t i = 0; i < u.size(); ++i) {
        plus[i] += h*v[i];
        minus[i] -= h*v[i];
    }

    const auto gp = problem.evaluate(plus).constraints;
    const auto gm = problem.evaluate(minus).constraints;

    require(gp.size() == analytic.size(),
            "Constraint finite-difference dimensions");

    double maximum_error = 0.0;

    for (std::size_t i = 0; i < gp.size(); ++i) {
        const double numerical =
            (gp[i]-gm[i])/(2.0*h);

        const double scale = std::max({
            1.0,
            std::abs(numerical),
            std::abs(analytic[i])
        });

        maximum_error = std::max(
            maximum_error,
            std::abs(numerical-analytic[i])/scale);
    }

    std::cout
        << "Maximum scaled constraint FD error: "
        << maximum_error << '\n';

    require(maximum_error < 1.0e-5,
            "Scaled constraint finite differences failed");
}

void test_all_fixed()
{
    auto cfg = configuration();

    for (auto& variable : cfg.variables)
        variable.free = false;

    const FullSpaceProblem problem(cfg);

    require(problem.free_variable_count() == 0,
            "All-fixed problem has free variables");

    const auto result = problem.evaluate({});

    require(result.objective_gradient.empty(),
            "All-fixed objective gradient not empty");

    require(!result.constraints.empty(),
            "All-fixed constraints missing");

    require(problem.constraint_jvp({}, {}).size() ==
            problem.constraint_count(),
            "All-fixed JVP dimension mismatch");

    for (double value : problem.constraint_jvp({}, {}))
        require(value == 0.0,
                "All-fixed JVP nonzero");

    require(problem.constraint_vjp(
                {}, dual(problem.constraint_count())).empty(),
            "All-fixed VJP not empty");

    require(problem.coordinates(
                problem.trajectory({})).empty(),
            "All-fixed coordinate recovery not empty");
}

void test_invalid_configuration()
{
    auto cfg = configuration();

    cfg.variables.pop_back();

    expect_invalid([&] {
        FullSpaceProblem problem(cfg);
    });

    cfg = configuration();
    cfg.variables[0].scale = 0.0;

    expect_invalid([&] {
        FullSpaceProblem problem(cfg);
    });

    cfg = configuration();
    cfg.variables[0].scale = -1.0;

    expect_invalid([&] {
        FullSpaceProblem problem(cfg);
    });

    cfg = configuration();
    cfg.variables[0].scale =
        std::numeric_limits<double>::infinity();

    expect_invalid([&] {
        FullSpaceProblem problem(cfg);
    });

    cfg = configuration();
    cfg.constraint_scales[0] = 0.0;

    expect_invalid([&] {
        FullSpaceProblem problem(cfg);
    });

    cfg = configuration();
    cfg.measurements.pop_back();

    expect_invalid([&] {
        FullSpaceProblem problem(cfg);
    });

    const FullSpaceProblem problem(configuration());

    expect_invalid([&] {
        problem.evaluate({1.0});
    });

    expect_invalid([&] {
        problem.constraint_vjp({}, {});
    });
}

} // namespace

int main()
{
    try {
        test_round_trip();
        test_reference_evaluation();

        for (int n : {2, 5, 10})
            for (int steps : {1, 2, 4})
                test_transpose(n, steps);

        test_objective_gradient();
        test_constraint_finite_difference();
        test_all_fixed();
        test_invalid_configuration();

        std::cout
            << "PASS E2.4b.10 full-space problem\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E2.4b.10: "
            << error.what() << '\n';

        return 1;
    }
}
