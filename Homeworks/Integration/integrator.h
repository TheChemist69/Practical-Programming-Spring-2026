// "integrator.h" header file.
// Adaptive quadrature: plain open 4-point (Part A) and Clenshaw-Curtis
// variable-transformation variants (Part B).
#pragma once

#include <cstddef>
#include <functional>
#include <string>

namespace pp {

// Outcome status for a full adaptive integration call.
enum class IntegratorStatus {
    success,
    max_depth_reached,
    interval_too_small,
    non_finite_evaluation
};

// User-tunable controls for the adaptive integrator.
struct AdaptiveOptions {
    double acc = 1e-6;          // absolute tolerance: accept when error < acc
    double eps = 1e-6;          // relative tolerance: accept when error < eps * |value|
    int max_depth = 100000;     // maximum recursion depth (safety guard)
    double min_interval = 0.0;  // stop recursing when interval width falls below this
};

// Result bundle returned by the adaptive integrator.
struct AdaptiveResult {
    double value = 0.0;                                       // estimated integral
    double estimated_error = 0.0;                             // accumulated error estimate
    std::size_t evaluations = 0;                              // total integrand calls
    IntegratorStatus status = IntegratorStatus::success;      // completion state
    int max_depth_reached = 0;                                // deepest recursion level seen
};

// Part A -- Integrates f on [a,b] using recursive open 4-point adaptive quadrature.
AdaptiveResult integrate_open4_adaptive(
        const std::function<double(double)>& f,
        double a,
        double b,
        const AdaptiveOptions& options = {});

// Part B -- Clenshaw-Curtis integrator on [a,b].
// Applies the substitution x = (a+b)/2 + (b-a)/2 * cos(theta), theta in [0,pi],
// which maps endpoint singularities (e.g. 1/sqrt(b-x)) to smooth sin(theta) factors.
// The transformed integrand is then passed to integrate_open4_adaptive.
AdaptiveResult integrate_clenshaw_curtis(
        const std::function<double(double)>& f,
        double a,
        double b,
        const AdaptiveOptions& options = {});

// Part B -- Semi-infinite integral on [a, +inf).
// Uses the substitution t = x/(1+x-a), x = a + t/(1-t), jacobian = 1/(1-t)^2,
// mapping [a,+inf) to [0,1), then applies integrate_open4_adaptive on (0,1).
AdaptiveResult integrate_semi_infinite(
        const std::function<double(double)>& f,
        double a,
        const AdaptiveOptions& options = {});

// Part B -- Doubly infinite integral on (-inf, +inf).
// Uses the substitution x = t/(1-t^2), jacobian = (1+t^2)/(1-t^2)^2,
// mapping (-inf,+inf) to (-1,1), then applies integrate_open4_adaptive on (-1,1).
AdaptiveResult integrate_doubly_infinite(
        const std::function<double(double)>& f,
        const AdaptiveOptions& options = {});

// Converts status enum values to short printable text labels.
const char* integrator_status_cstr(IntegratorStatus status);

} // namespace pp
