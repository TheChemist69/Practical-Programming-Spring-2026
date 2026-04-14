// "main.cc" implementation file.
// Orchestration for recursive adaptive integration homework tasks.

#include "benchmarks.h"
#include "erf_integral.h"
#include "reporting.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace pp {

namespace {

// Runs the assignment benchmark integrals and prints pass/fail diagnostics.
void run_benchmark_section() {
    std::puts("=== Recursive open 4-point adaptive integrator ===");
    std::puts("Benchmark integrals on [0,1] with mixed absolute/relative tolerance");

    AdaptiveOptions options;
    options.acc = 1e-6;
    options.eps = 1e-6;
    options.max_depth = 100000;
    options.min_interval = 0.0;

    const std::vector<BenchmarkResult> rows = run_benchmarks(options);
    write_benchmark_results("benchmark_results.data", rows);

    int passes = 0;
    for (const BenchmarkResult& row : rows) {
        if (row.passed) {
            passes++;
        }
        std::printf(
                "  %-28s Q=% .12g exact=% .12g |err|=%.3e est=%.3e tol=%.3e calls=%zu status=%s %s\n",
                row.name.c_str(),
                row.result.value,
                row.exact,
                row.abs_error,
                row.result.estimated_error,
                row.tolerance_goal,
                row.result.evaluations,
                integrator_status_cstr(row.result.status),
                row.passed ? "PASS" : "FAIL");
    }

    std::printf("Summary: %d/%zu benchmark checks passed\n", passes, rows.size());
    std::puts("Note: int_0_1_sqrt_1_minus_x2 exact value is pi/4 (checked analytically).\n");
}

// Computes erf data and prints the requested erf(1) comparison.
void run_erf_section() {
    std::puts("=== Error function via integral representation ===");

    AdaptiveOptions curve_options;
    curve_options.acc = 1e-8;
    curve_options.eps = 1e-8;
    curve_options.max_depth = 100000;
    curve_options.min_interval = 0.0;

    const std::vector<ErfCurvePoint> curve = sample_erf_curve(-3.0, 3.0, 0.05, curve_options);
    const std::vector<std::pair<double, double>> tabulated = erf_tabulated_values();

    write_erf_curve_data("erf_curve.data", curve);
    write_erf_tabulated_data("erf_tabulated.data", tabulated);

    const ErfEvaluation erf1 = erf_via_integral(1.0, curve_options);
    const double reference = erf1_reference_value();
    const double abs_error = std::abs(erf1.value - reference);

    std::printf("erf(1) numeric  = %.17g\n", erf1.value);
    std::printf("erf(1) reference= %.17g\n", reference);
    std::printf("|difference|    = %.3e\n", abs_error);
    std::printf("estimated error = %.3e\n", erf1.estimated_error);
    std::printf("function calls  = %zu\n", erf1.evaluations);
    std::printf("status          = %s\n\n", integrator_status_cstr(erf1.status));
}

// Builds the decreasing-acc convergence table for erf(1).
void run_erf_sweep_section() {
    std::puts("=== erf(1) convergence with eps=0 and decreasing acc ===");

    AdaptiveOptions sweep_options;
    sweep_options.max_depth = 100000;
    sweep_options.min_interval = 0.0;

    const std::vector<ErfSweepPoint> rows = erf1_accuracy_sweep(1e-1, 10, sweep_options);
    write_erf_sweep_data("erf_acc_sweep.data", rows);

    std::puts("  acc          |erf(1)-exact|    est_error        calls   status");
    for (const ErfSweepPoint& row : rows) {
        std::printf("  %-10.1e %.3e        %.3e   %8zu   %s\n",
                row.acc,
                row.abs_error,
                row.estimated_error,
                row.evaluations,
                integrator_status_cstr(row.status));
    }
    std::puts("");
}

} // namespace

// Executes sections in assignment order and lists generated outputs.
int run() {
    run_benchmark_section();
    run_erf_section();
    run_erf_sweep_section();

    std::puts("Data files written:");
    std::puts("  benchmark_results.data");
    std::puts("  erf_curve.data");
    std::puts("  erf_tabulated.data");
    std::puts("  erf_acc_sweep.data");

    return EXIT_SUCCESS;
}

} // namespace pp

// Program entry point.
int main() {
    return pp::run();
}
