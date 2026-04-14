// "erf_integral.cc" implementation file.
// Error function via adaptive integral evaluation.

#include "erf_integral.h"

#include <cmath>
#include <vector>

namespace pp {

namespace {

// External reference used for convergence and plotting checks.
constexpr double kErf1Reference = 0.84270079294971486934;

// Converts an integral estimate to erf(z) and propagates error scaling.
ErfEvaluation map_integral_to_erf(double z, const AdaptiveResult& integral_result, bool transformed_branch) {
    const double pi = std::acos(-1.0);
    const double factor = 2.0 / std::sqrt(pi);

    ErfEvaluation out;
    out.z = z;
    out.evaluations = integral_result.evaluations;
    out.status = integral_result.status;

    if (transformed_branch) {
        out.value = 1.0 - factor * integral_result.value;
    } else {
        out.value = factor * integral_result.value;
    }
    out.estimated_error = factor * integral_result.estimated_error;
    return out;
}

} // namespace

// Returns assignment reference value for erf(1).
double erf1_reference_value() {
    return kErf1Reference;
}

// Piecewise erf(z): odd symmetry, direct integral, or transformed tail integral.
ErfEvaluation erf_via_integral(double z, const AdaptiveOptions& options) {
    // Use odd symmetry to map negative z to positive side.
    if (z < 0.0) {
        ErfEvaluation pos = erf_via_integral(-z, options);
        pos.z = z;
        pos.value = -pos.value;
        return pos;
    }

    // Fast path at the origin.
    if (z == 0.0) {
        return ErfEvaluation{};
    }

    // Direct representation on [0,1].
    if (z <= 1.0) {
        AdaptiveResult res = integrate_open4_adaptive(
                [](double x) { return std::exp(-x * x); },
                0.0,
                z,
                options);
        return map_integral_to_erf(z, res, false);
    }

    // Tail representation on t in (0,1] for z > 1.
    AdaptiveResult res = integrate_open4_adaptive(
            [z](double t) {
                const double inv_t = 1.0 / t;
                const double u = z + (1.0 - t) * inv_t;
                return std::exp(-u * u) * inv_t * inv_t;
            },
            0.0,
            1.0,
            options);
    return map_integral_to_erf(z, res, true);
}

        // Uniformly samples erf(z) and stores error against std::erf.
std::vector<ErfCurvePoint> sample_erf_curve(
        double z_min,
        double z_max,
        double dz,
        const AdaptiveOptions& options) {
    std::vector<ErfCurvePoint> rows;

    if (dz <= 0.0 || z_max < z_min) {
        return rows;
    }

    const double guard = 0.5 * dz;
    for (double z = z_min; z <= z_max + guard; z += dz) {
        ErfEvaluation eval = erf_via_integral(z, options);

        ErfCurvePoint row;
        row.z = z;
        row.value = eval.value;
        row.reference = std::erf(z);
        row.abs_error = std::abs(row.value - row.reference);
        row.estimated_error = eval.estimated_error;
        row.evaluations = eval.evaluations;
        row.status = eval.status;
        rows.push_back(row);
    }

    return rows;
}

// Generates erf(1) accuracy rows for acc, acc/10, acc/100, ...
std::vector<ErfSweepPoint> erf1_accuracy_sweep(
        double acc_start,
        int levels,
        const AdaptiveOptions& base_options) {
    std::vector<ErfSweepPoint> rows;
    if (acc_start <= 0.0 || levels <= 0) {
        return rows;
    }

    rows.reserve(static_cast<std::size_t>(levels));

    for (int i = 0; i < levels; i++) {
        const double acc = acc_start * std::pow(0.1, i);

        AdaptiveOptions opts = base_options;
        opts.acc = acc;
        opts.eps = 0.0;

        ErfEvaluation eval = erf_via_integral(1.0, opts);

        ErfSweepPoint row;
        row.acc = acc;
        row.value = eval.value;
        row.abs_error = std::abs(eval.value - kErf1Reference);
        row.estimated_error = eval.estimated_error;
        row.evaluations = eval.evaluations;
        row.status = eval.status;
        rows.push_back(row);
    }

    return rows;
}

// Tabulated values used in assignment-style comparison plots.
std::vector<std::pair<double, double>> erf_tabulated_values() {
    return {
        {0.0, 0.0},
        {0.1, 0.112462916},
        {0.2, 0.222702589},
        {0.3, 0.328626759},
        {0.4, 0.428392355},
        {0.5, 0.520499878},
        {0.6, 0.603856091},
        {0.7, 0.677801194},
        {0.8, 0.742100965},
        {0.9, 0.796908212},
        {1.0, 0.842700793},
        {1.5, 0.966105146},
        {2.0, 0.995322265},
        {2.5, 0.999593048},
        {3.0, 0.999977910},
    };
}

} // namespace pp
