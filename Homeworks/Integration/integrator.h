// "integrator.h" header file.
// Recursive open 4-point adaptive quadrature declarations.
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

// User-tunable controls for tolerance and recursion safety.
struct AdaptiveOptions {
    double acc = 1e-6;
    double eps = 1e-6;
    int max_depth = 100000;
    double min_interval = 0.0;
};

// Result bundle returned by the adaptive integrator.
struct AdaptiveResult {
    double value = 0.0;
    double estimated_error = 0.0;
    std::size_t evaluations = 0;
    IntegratorStatus status = IntegratorStatus::success;
    int max_depth_reached = 0;
};

// Integrates f on [a,b] using recursive open 4-point adaptive quadrature.
AdaptiveResult integrate_open4_adaptive(
        const std::function<double(double)>& f,
        double a,
        double b,
        const AdaptiveOptions& options = {});

// Converts status enum values to short printable text labels.
const char* integrator_status_cstr(IntegratorStatus status);

} // namespace pp
