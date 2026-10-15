#include "tpsf/dense_full_space_subproblem.hpp"
#include "tpsf/dense_constrained_least_squares.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

using namespace tpsf;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

double relative_error(
    const Eigen::VectorXd& a,
    const Eigen::VectorXd& b)
{
    return (a-b).norm() /
        std::max({1.0, a.norm(), b.norm()});
}

FullSpaceProblemConfiguration make_configuration(
    bool all_fixed = false,
    bool no_measurements = false)
{
    FullSpaceProblemConfiguration config;

    constexpr int n = 2;
    constexpr std::size_t levels = 3;
    constexpr std::size_t width = 2*n+2;
    constexpr std::size_t block = width+6;

    config.reference.n = n;

    for (std::size_t k = 0; k < levels; ++k) {
        std::vector<double> x(width);

        for (std::size_t i = 0; i < width; ++i) {
            const bool pressure =
                i == 0 || (i >= 3 && i % 2 == 1);

            x[i] = pressure
                ? 5.0e6-1000.0*i+200.0*k
                : 10.0+0.1*i+0.05*k;
        }

        config.reference.states.push_back(x);

        PhysicalBoundaryInputs q;
        q.inlet_pressure = 5.0e6+200.0*k;
        q.outlet_flow = 10.0+0.1*k;

        config.reference.inputs.push_back(q);

        BoundaryVector b{
            q.inlet_pressure+100.0,
            x[1]+0.1,
            x[2*n+1]+100.0,
            q.outlet_flow+0.1
        };

        config.reference.reconstructed.push_back(b);

        BoundaryMeasurement measurement;

        for (std::size_t c = 0; c < 4; ++c) {
            const bool pressure = c == 0 || c == 2;

            measurement.available[c] =
                !no_measurements && (k+c)%3 != 0;

            measurement.sigma[c] =
                pressure ? 200.0 : 0.2;

            measurement.observed[c] =
                b[c]-0.5*measurement.sigma[c];
        }

        config.measurements.push_back(measurement);
    }

    config.parameters = pipe_sim::TransientParameters{
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
    };

    config.variables.resize(levels*block);

    for (std::size_t k = 0; k < levels; ++k) {
        for (std::size_t i = 0; i < block; ++i) {
            const auto index = k*block+i;

            const bool pressure =
                i == 0 ||
                (i >= 3 && i < width && i%2 == 1) ||
                i == width ||
                i == width+2 ||
                i == width+4;

            config.variables[index].scale =
                pressure ? 1.0e5 : 5.0;

            config.variables[index].free = !all_fixed;

            // Prescribed inlet pressure remains fixed.
            if (i == width)
                config.variables[index].free = false;
        }
    }

    const std::size_t constraints =
        (levels-1)*width+4*levels;

    config.constraint_scales.resize(constraints);

    for (std::size_t i = 0; i < constraints; ++i)
        config.constraint_scales[i] =
            i%3 == 0 ? 1.0e-5 :
            i%3 == 1 ? 0.1 : 2.0;

    return config;
}

void verify_assembly(
    const FullSpaceProblem& problem)
{
    Eigen::VectorXd u = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(
            problem.free_variable_count()));

    for (Eigen::Index j = 0; j < u.size(); ++j)
        u[j] = 0.01*std::sin(0.3*static_cast<double>(j));

    const auto assembled =
        assemble_dense_full_space_subproblem(problem, u);

    const auto coordinates =
        std::vector<double>(u.data(), u.data()+u.size());

    const auto evaluation = problem.evaluate(coordinates);

    require(
        assembled.constraint_jacobian.rows() ==
            static_cast<Eigen::Index>(
                problem.constraint_count()),
        "Constraint row count");

    require(
        assembled.constraint_jacobian.cols() == u.size(),
        "Constraint column count");

    require(
        assembled.measurement_jacobian.cols() == u.size(),
        "Measurement column count");

    const double objective =
        0.5*assembled.measurement_residual.squaredNorm();

    require(
        std::abs(objective-evaluation.measurement.value)
            <= 1e-11*std::max(1.0, objective),
        "Measurement objective identity");

    const Eigen::VectorXd gradient =
        assembled.measurement_jacobian.transpose() *
        assembled.measurement_residual;

    Eigen::VectorXd expected_gradient(
        evaluation.objective_gradient.size());

    for (std::size_t i = 0;
         i < evaluation.objective_gradient.size(); ++i)
        expected_gradient[static_cast<Eigen::Index>(i)] =
            evaluation.objective_gradient[i];

    require(
        relative_error(gradient, expected_gradient) < 1e-11,
        "Measurement gradient identity");

    Eigen::VectorXd expected_constraints(
        evaluation.constraints.size());

    for (std::size_t i = 0;
         i < evaluation.constraints.size(); ++i)
        expected_constraints[static_cast<Eigen::Index>(i)] =
            evaluation.constraints[i];

    require(
        relative_error(
            assembled.constraint_residual,
            expected_constraints) < 1e-12,
        "Constraint residual identity");

    Eigen::VectorXd dual =
        Eigen::VectorXd::Zero(
            assembled.constraint_jacobian.rows());

    for (Eigen::Index i = 0; i < dual.size(); ++i)
        dual[i] = std::cos(0.2*static_cast<double>(i));

    const auto dual_std =
        std::vector<double>(
            dual.data(), dual.data()+dual.size());

    const auto vjp =
        problem.constraint_vjp(coordinates, dual_std);

    Eigen::VectorXd expected_vjp(vjp.size());

    for (std::size_t i = 0; i < vjp.size(); ++i)
        expected_vjp[static_cast<Eigen::Index>(i)] = vjp[i];

    require(
        relative_error(
            assembled.constraint_jacobian.transpose()*dual,
            expected_vjp) < 1e-10,
        "Constraint VJP transpose identity");

    // Independently check selected columns against centered
    // finite differences of the actual nonlinear constraints.
    if (u.size() > 0) {
        constexpr double h = 1e-5;

        const Eigen::Index selected[] = {
            0, u.size()/2, u.size()-1
        };

        for (const auto j : selected) {
            auto plus = u;
            auto minus = u;

            plus[j] += h;
            minus[j] -= h;

            const auto plus_eval = problem.evaluate(
                std::vector<double>(
                    plus.data(), plus.data()+plus.size()));

            const auto minus_eval = problem.evaluate(
                std::vector<double>(
                    minus.data(), minus.data()+minus.size()));

            Eigen::VectorXd finite_difference(
                plus_eval.constraints.size());

            for (std::size_t i = 0;
                 i < plus_eval.constraints.size(); ++i) {
                finite_difference[
                    static_cast<Eigen::Index>(i)] =
                    (plus_eval.constraints[i] -
                     minus_eval.constraints[i])/(2*h);
            }

            require(
                relative_error(
                    finite_difference,
                    assembled.constraint_jacobian.col(j))
                    < 1e-5,
                "Constraint finite-difference check");

            const auto plus_residual =
                problem.measurement_residual(
                    std::vector<double>(
                        plus.data(), plus.data()+plus.size()));

            const auto minus_residual =
                problem.measurement_residual(
                    std::vector<double>(
                        minus.data(), minus.data()+minus.size()));

            Eigen::VectorXd measurement_fd(
                plus_residual.size());

            for (std::size_t i = 0;
                 i < plus_residual.size(); ++i) {
                measurement_fd[
                    static_cast<Eigen::Index>(i)] =
                    (plus_residual[i] -
                     minus_residual[i])/(2*h);
            }

            require(
                relative_error(
                    measurement_fd,
                    assembled.measurement_jacobian.col(j))
                    < 1e-5,
                "Measurement finite-difference check");
        }
    }
}

} // namespace

int main()
{
    try {
        const FullSpaceProblem ordinary(
            make_configuration());

        verify_assembly(ordinary);

        const FullSpaceProblem no_measurements(
            make_configuration(false, true));

        verify_assembly(no_measurements);

        const auto empty =
            assemble_dense_full_space_subproblem(
                no_measurements,
                Eigen::VectorXd::Zero(
                    no_measurements.free_variable_count()));

        require(
            empty.measurement_jacobian.rows() == 0,
            "Missing measurements created Jacobian rows");

        const FullSpaceProblem all_fixed(
            make_configuration(true, false));

        verify_assembly(all_fixed);

        const auto fixed =
            assemble_dense_full_space_subproblem(
                all_fixed, Eigen::VectorXd(0));

        require(
            fixed.constraint_jacobian.cols() == 0,
            "Fixed variables created Jacobian columns");

        require(
            fixed.measurement_jacobian.cols() == 0,
            "Fixed variables created measurement columns");

        std::cout
            << "PASS E2.4b.11b.1 dense full-space assembly\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL E2.4b.11b.1: "
            << error.what() << '\n';

        return 1;
    }
}
