// "erf_integral.h" header file.
// Error function via adaptive-integral representation.
#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "integrator.h"

namespace pp {

// Single erf(z) evaluation with integration diagnostics.
struct ErfEvaluation {
    double z = 0.0;
    double value = 0.0;
    double estimated_error = 0.0;
    std::size_t evaluations = 0;
    IntegratorStatus status = IntegratorStatus::success;
};

// One point on an erf curve including reference comparison.
struct ErfCurvePoint {
    double z = 0.0;
    double value = 0.0;
    double reference = 0.0;
    double abs_error = 0.0;
    double estimated_error = 0.0;
    std::size_t evaluations = 0;
    IntegratorStatus status = IntegratorStatus::success;
};

// One row in the erf(1) accuracy-vs-acc sweep.
struct ErfSweepPoint {
    double acc = 0.0;
    double value = 0.0;
    double abs_error = 0.0;
    double estimated_error = 0.0;
    std::size_t evaluations = 0;
    IntegratorStatus status = IntegratorStatus::success;
};

    // Computes erf(z) from the assignment's piecewise integral formula.
ErfEvaluation erf_via_integral(double z, const AdaptiveOptions& options = {});

    // Samples erf(z) on a uniform grid and compares to std::erf.
std::vector<ErfCurvePoint> sample_erf_curve(
        double z_min,
        double z_max,
        double dz,
        const AdaptiveOptions& options = {});

    // Sweeps acc values for erf(1) with eps forced to zero.
std::vector<ErfSweepPoint> erf1_accuracy_sweep(
        double acc_start,
        int levels,
        const AdaptiveOptions& base_options = {});

    // Tabulated points used for visual comparison in plots.
std::vector<std::pair<double, double>> erf_tabulated_values();

    // High-precision reference value for erf(1).
double erf1_reference_value();

} // namespace pp
