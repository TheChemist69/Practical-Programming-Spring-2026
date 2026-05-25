// "integrator.cc" implementation file.
// Recursive open 4-point adaptive quadrature (Part A) and
// Clenshaw-Curtis / infinite-limit transformations (Part B).

#include "integrator.h"

#include <cmath>
#include <limits>

namespace pp {

namespace {

// Internal return type used while descending the recursion tree.
struct RecursiveResult {
    double value = 0.0;
    double error = 0.0;
    std::size_t evaluations = 0;
    IntegratorStatus status = IntegratorStatus::success;
    int max_depth_reached = 0;
};

// Prioritize the most severe status seen in child branches.
IntegratorStatus merge_status(IntegratorStatus a, IntegratorStatus b) {
    if (a == IntegratorStatus::non_finite_evaluation || b == IntegratorStatus::non_finite_evaluation) {
        return IntegratorStatus::non_finite_evaluation;
    }
    if (a == IntegratorStatus::max_depth_reached || b == IntegratorStatus::max_depth_reached) {
        return IntegratorStatus::max_depth_reached;
    }
    if (a == IntegratorStatus::interval_too_small || b == IntegratorStatus::interval_too_small) {
        return IntegratorStatus::interval_too_small;
    }
    return IntegratorStatus::success;
}

// Recursive core implementing open 4-point adaptive subdivision.
RecursiveResult integrate_recursive(
        const std::function<double(double)>& f,
        double a,
        double b,
        double acc,
        double eps,
        double f2,
        double f3,
        int depth,
        const AdaptiveOptions& options) {
    RecursiveResult out;
    out.max_depth_reached = depth;

    // Segment size and reusable inner nodes for the open rule.
    const double h = b - a;
    const double x2 = a + 2.0 * h / 6.0;
    const double x3 = a + 4.0 * h / 6.0;

    // Reuse parent evaluations when available.
    if (std::isnan(f2)) {
        f2 = f(x2);
        out.evaluations++;
    }
    if (std::isnan(f3)) {
        f3 = f(x3);
        out.evaluations++;
    }

    // Evaluate the two outer open points for this interval.
    const double x1 = a + h / 6.0;
    const double x4 = a + 5.0 * h / 6.0;
    const double f1 = f(x1);
    const double f4 = f(x4);
    out.evaluations += 2;

    if (!std::isfinite(f1) || !std::isfinite(f2) || !std::isfinite(f3) || !std::isfinite(f4)) {
        out.value = std::numeric_limits<double>::quiet_NaN();
        out.error = std::numeric_limits<double>::infinity();
        out.status = IntegratorStatus::non_finite_evaluation;
        return out;
    }

    // Embedded pair: Q is higher-order, q is lower-order.
    const double q_high = (2.0 * f1 + f2 + f3 + 2.0 * f4) / 6.0 * h;
    const double q_low = (f1 + f2 + f3 + f4) / 4.0 * h;
    const double err = std::abs(q_high - q_low);
    const double tol = acc + eps * std::abs(q_high);

    out.value = q_high;
    out.error = err;

    if (!std::isfinite(q_high) || !std::isfinite(err)) {
        out.value = std::numeric_limits<double>::quiet_NaN();
        out.error = std::numeric_limits<double>::infinity();
        out.status = IntegratorStatus::non_finite_evaluation;
        return out;
    }

    // Accept this segment if local error meets local tolerance.
    if (err <= tol) {
        return out;
    }

    // Stop when recursion limits or floating-point limits are reached.
    if (depth >= options.max_depth) {
        out.status = IntegratorStatus::max_depth_reached;
        return out;
    }

    const double mid = 0.5 * (a + b);
    if (mid == a || mid == b || (options.min_interval > 0.0 && std::abs(h) <= options.min_interval)) {
        out.status = IntegratorStatus::interval_too_small;
        return out;
    }

    // Split and distribute absolute tolerance between child intervals.
    const double child_acc = acc / std::sqrt(2.0);

    RecursiveResult left = integrate_recursive(f, a, mid, child_acc, eps, f1, f2, depth + 1, options);
    RecursiveResult right = integrate_recursive(f, mid, b, child_acc, eps, f3, f4, depth + 1, options);

    // Combine child integrals and child error estimates.
    out.value = left.value + right.value;
    out.error = std::sqrt(left.error * left.error + right.error * right.error);
    out.evaluations += left.evaluations + right.evaluations;
    out.status = merge_status(left.status, right.status);
    out.max_depth_reached = (left.max_depth_reached > right.max_depth_reached)
            ? left.max_depth_reached
            : right.max_depth_reached;

    return out;
}

} // namespace

// Public wrapper that validates inputs and starts recursion.
AdaptiveResult integrate_open4_adaptive(
        const std::function<double(double)>& f,
        double a,
        double b,
        const AdaptiveOptions& options) {
    AdaptiveResult out;

    if (a == b) {
        return out;
    }

    if (options.max_depth < 0 || options.acc < 0.0 || options.eps < 0.0 || options.min_interval < 0.0) {
        out.value = std::numeric_limits<double>::quiet_NaN();
        out.estimated_error = std::numeric_limits<double>::infinity();
        out.status = IntegratorStatus::non_finite_evaluation;
        return out;
    }

    const double missing = std::numeric_limits<double>::quiet_NaN();
    RecursiveResult rec = integrate_recursive(
            f,
            a,
            b,
            options.acc,
            options.eps,
            missing,
            missing,
            0,
            options);

    out.value = rec.value;
    out.estimated_error = rec.error;
    out.evaluations = rec.evaluations;
    out.status = rec.status;
    out.max_depth_reached = rec.max_depth_reached;
    return out;
}

// Part B -- Clenshaw-Curtis integrator on [a,b].
// Substitution: x = m + r*cos(theta),  theta in [0, pi]
// where m = (a+b)/2, r = (b-a)/2.
// dx = -r*sin(theta) dtheta, limits flip:
//   integral_a^b f(x) dx = integral_0^pi f(m + r*cos(t)) * r*sin(t) dt
// This is smooth even when f has integrable endpoint singularities such as
// 1/sqrt(b-x) or log(x-a), because sin(theta)->0 cancels the divergence.
AdaptiveResult integrate_clenshaw_curtis(
        const std::function<double(double)>& f,
        double a,
        double b,
        const AdaptiveOptions& options) {
    const double m = 0.5 * (a + b);
    const double r = 0.5 * (b - a);
    // Transformed integrand in the theta variable.
    const auto g = [&](double theta) -> double {
        const double x = m + r * std::cos(theta);
        return f(x) * std::sin(theta) * r;
    };
    return integrate_open4_adaptive(g, 0.0, M_PI, options);
}

// Part B -- Semi-infinite integrator on [a, +inf).
// Substitution: t = (x-a)/(1+(x-a)),  x = a + t/(1-t),  t in [0,1)
// Jacobian: dx/dt = 1/(1-t)^2.
// The transformed integrand is integrated on (0,1).
AdaptiveResult integrate_semi_infinite(
        const std::function<double(double)>& f,
        double a,
        const AdaptiveOptions& options) {
    const auto g = [&](double t) -> double {
        // Guard against t reaching exactly 1 (open upper end).
        if (t >= 1.0) return 0.0;
        const double one_minus_t = 1.0 - t;
        const double x = a + t / one_minus_t;
        const double jacobian = 1.0 / (one_minus_t * one_minus_t);
        return f(x) * jacobian;
    };
    return integrate_open4_adaptive(g, 0.0, 1.0, options);
}

// Part B -- Doubly infinite integrator on (-inf, +inf).
// Substitution: x = t/(1-t^2),  t in (-1,1)
// Jacobian: dx/dt = (1+t^2)/(1-t^2)^2.
// The transformed integrand is integrated on (-1,1).
AdaptiveResult integrate_doubly_infinite(
        const std::function<double(double)>& f,
        const AdaptiveOptions& options) {
    const auto g = [&](double t) -> double {
        // Guard against t reaching ±1 (open endpoints).
        if (std::abs(t) >= 1.0) return 0.0;
        const double t2 = t * t;
        const double denom = 1.0 - t2;
        const double x = t / denom;
        const double jacobian = (1.0 + t2) / (denom * denom);
        return f(x) * jacobian;
    };
    return integrate_open4_adaptive(g, -1.0, 1.0, options);
}

// Human-readable labels used in logs and output files.
const char* integrator_status_cstr(IntegratorStatus status) {
    switch (status) {
        case IntegratorStatus::success:
            return "success";
        case IntegratorStatus::max_depth_reached:
            return "max_depth_reached";
        case IntegratorStatus::interval_too_small:
            return "interval_too_small";
        case IntegratorStatus::non_finite_evaluation:
            return "non_finite_evaluation";
    }
    return "unknown";
}

} // namespace pp
