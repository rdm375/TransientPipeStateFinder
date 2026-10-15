#include "tpsf/minres.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace tpsf {
namespace {

double dot(
    const KrylovVector& a,
    const KrylovVector& b)
{
    if (a.size() != b.size())
        throw std::invalid_argument("Krylov vector size mismatch");

    double result = 0.0;

    for (std::size_t i = 0; i < a.size(); ++i)
        result += a[i]*b[i];

    return result;
}

double norm(const KrylovVector& v)
{
    // Avoid overflow/underflow in sum of squares.
    double scale = 0.0;
    double ssq = 1.0;

    for (double value : v) {
        if (!std::isfinite(value))
            throw std::invalid_argument("Nonfinite Krylov vector");

        const double a = std::abs(value);

        if (a == 0.0)
            continue;

        if (scale < a) {
            ssq = 1.0 + ssq*(scale/a)*(scale/a);
            scale = a;
        }
        else {
            ssq += (a/scale)*(a/scale);
        }
    }

    return scale == 0.0 ? 0.0 : scale*std::sqrt(ssq);
}

KrylovVector checked_apply(
    const SymmetricOperator& apply,
    const KrylovVector& x)
{
    auto y = apply(x);

    if (y.size() != x.size())
        throw std::invalid_argument(
            "Krylov operator returned incorrect dimension");

    (void)norm(y);
    return y;
}

double true_residual(
    const SymmetricOperator& apply,
    const KrylovVector& x,
    const KrylovVector& rhs)
{
    auto ax = checked_apply(apply, x);

    for (std::size_t i = 0; i < ax.size(); ++i)
        ax[i] = rhs[i]-ax[i];

    return norm(ax);
}

} // namespace

MinresResult solve_minres(
    const SymmetricOperator& apply,
    const KrylovVector& rhs,
    const MinresOptions& options)
{
    if (!apply)
        throw std::invalid_argument("Missing Krylov operator");

    if (!std::isfinite(options.relative_tolerance) ||
        !std::isfinite(options.absolute_tolerance) ||
        options.relative_tolerance < 0.0 ||
        options.absolute_tolerance < 0.0 ||
        options.max_iterations == 0) {
        throw std::invalid_argument("Invalid MINRES options");
    }

    const std::size_t n = rhs.size();

    MinresResult result;
    result.solution.assign(n, 0.0);

    const double beta1 = norm(rhs);

    result.initial_residual_norm = beta1;
    result.final_residual_norm = beta1;

    if (beta1 == 0.0) {
        result.status = MinresStatus::converged;
        return result;
    }

    const double target = std::max(
        options.absolute_tolerance,
        options.relative_tolerance*beta1);

    KrylovVector r1(rhs);
    KrylovVector r2(rhs);

    KrylovVector w(n, 0.0);
    KrylovVector w1(n, 0.0);
    KrylovVector w2(n, 0.0);

    double beta = beta1;
    double oldb = 0.0;

    double dbar = 0.0;
    double epsln = 0.0;
    double phibar = beta1;

    double cs = -1.0;
    double sn = 0.0;

    for (std::size_t iteration = 1;
         iteration <= options.max_iterations;
         ++iteration) {

        if (!(beta > 0.0) || !std::isfinite(beta)) {
            result.status = MinresStatus::breakdown;
            return result;
        }

        KrylovVector v(n);

        for (std::size_t i = 0; i < n; ++i)
            v[i] = r2[i]/beta;

        KrylovVector y = checked_apply(apply, v);

        if (iteration >= 2) {
            for (std::size_t i = 0; i < n; ++i)
                y[i] -= (beta/oldb)*r1[i];
        }

        const double alpha = dot(v, y);

        for (std::size_t i = 0; i < n; ++i)
            y[i] -= (alpha/beta)*r2[i];

        r1 = std::move(r2);
        r2 = std::move(y);

        oldb = beta;
        beta = norm(r2);

        const double oldeps = epsln;

        const double delta = cs*dbar + sn*alpha;
        const double gbar = sn*dbar - cs*alpha;

        epsln = sn*beta;
        dbar = -cs*beta;

        const double gamma = std::hypot(gbar, beta);

        if (!(gamma > 0.0) || !std::isfinite(gamma)) {
            result.status = MinresStatus::breakdown;
            result.iterations = iteration-1;
            return result;
        }

        cs = gbar/gamma;
        sn = beta/gamma;

        const double phi = cs*phibar;
        phibar = sn*phibar;

        w1 = std::move(w2);
        w2 = std::move(w);
        w.resize(n);

        for (std::size_t i = 0; i < n; ++i) {
            w[i] = (
                v[i] - oldeps*w1[i] - delta*w2[i]
            )/gamma;

            result.solution[i] += phi*w[i];
        }

        result.iterations = iteration;

        result.final_residual_norm = true_residual(
            apply, result.solution, rhs);

        if (!std::isfinite(result.final_residual_norm)) {
            result.status = MinresStatus::breakdown;
            return result;
        }

        if (result.final_residual_norm <= target) {
            result.status = MinresStatus::converged;
            return result;
        }

        // A completed Lanczos basis with nonzero residual
        // indicates an incompatible or numerically unresolved
        // system. Do not report convergence.
        if (beta == 0.0) {
            result.status = MinresStatus::breakdown;
            return result;
        }
    }

    result.status = MinresStatus::iteration_limit;
    return result;
}

} // namespace tpsf
