// Nikolaos Antonios Iakynthos Nemet
// This project was created using a combination of Claude Code and GitHub Copilot.

// "main.cc" implementation file.
// Orchestration for recursive adaptive integration homework tasks.
// Part A: plain open 4-point adaptive integrator and erf benchmark.
// Part B: Clenshaw-Curtis variable transformation and infinite-limit integrals.
// Part C: error estimate quality study (sweep over acc for erf(1)).

#include "benchmarks.h"
#include "erf_integral.h"
#include "integrator.h"
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
    run_benchmark_section();    // Part A: plain integrator benchmarks
    run_erf_section();          // Part A: erf via integral
    run_erf_sweep_section();    // Part C: error estimate quality study
    run_part_b_section();       // Part B: Clenshaw-Curtis + infinite limits

    std::puts("Data files written:");
    std::puts("  benchmark_results.data");
    std::puts("  erf_curve.data");
    std::puts("  erf_tabulated.data");
    std::puts("  erf_acc_sweep.data");
    std::puts("  cc_singular_comparison.data  (Part B: CC vs plain on singular integrand)");
    std::puts("  infinite_limits.data          (Part B: semi- and doubly-infinite results)");

    return EXIT_SUCCESS;
}

// Part B: demonstrates Clenshaw-Curtis and infinite-limit integration.
//
// Three tests are shown:
// 1. Singular integrand 1/sqrt(1-x^2) on (-1,1): exact = pi.
//    The plain integrator struggles because 1/sqrt(1-x^2) diverges at both
//    endpoints; the CC transformation smooths it to a constant pi on [0,pi].
// 2. Semi-infinite integral of exp(-x) on [0,+inf): exact = 1.
// 3. Doubly infinite integral of exp(-x^2) on (-inf,+inf): exact = sqrt(pi).
//
// For each case the number of function evaluations is compared between
// the plain integrator (on a truncated domain) and the transformed integrator.
void run_part_b_section() {
    std::puts("=== Part B: Clenshaw-Curtis and infinite-limit integrals ===");

    AdaptiveOptions opts;
    opts.acc = 1e-6;
    opts.eps = 1e-6;
    opts.max_depth = 100000;

    // ---- Test 1: integral_{-1}^{1} 1/sqrt(1-x^2) dx = pi ---------------
    // Plain integrator on (-1+eps, 1-eps) to avoid the endpoint divergences.
    const double a_sing = -1.0 + 1e-6, b_sing = 1.0 - 1e-6;
    const auto f_sing = [](double x) { return 1.0 / std::sqrt(1.0 - x*x); };

    AdaptiveResult plain_sing = integrate_open4_adaptive(f_sing, a_sing, b_sing, opts);
    AdaptiveResult cc_sing    = integrate_clenshaw_curtis(f_sing, -1.0, 1.0, opts);

    const double exact_sing = M_PI;
    std::printf("  Singular 1/sqrt(1-x^2) on (-1,1)  exact=pi=%.10g\n", exact_sing);
    std::printf("    plain (truncated):  value=%.10g  |err|=%.3e  evals=%zu\n",
                plain_sing.value, std::abs(plain_sing.value - exact_sing), plain_sing.evaluations);
    std::printf("    Clenshaw-Curtis:    value=%.10g  |err|=%.3e  evals=%zu\n",
                cc_sing.value, std::abs(cc_sing.value - exact_sing), cc_sing.evaluations);

    // Write comparison data for plotting.
    {
        FILE* fp = std::fopen("cc_singular_comparison.data", "w");
        if (fp) {
            std::fprintf(fp, "# method value abs_err evals\n");
            std::fprintf(fp, "plain %.15g %.3e %zu\n",
                         plain_sing.value, std::abs(plain_sing.value - exact_sing), plain_sing.evaluations);
            std::fprintf(fp, "cc    %.15g %.3e %zu\n",
                         cc_sing.value, std::abs(cc_sing.value - exact_sing), cc_sing.evaluations);
            std::fclose(fp);
        }
    }

    // ---- Test 2: integral_{0}^{+inf} exp(-x) dx = 1 --------------------
    const auto f_semi = [](double x) { return std::exp(-x); };
    AdaptiveResult semi = integrate_semi_infinite(f_semi, 0.0, opts);
    const double exact_semi = 1.0;
    std::printf("\n  Semi-infinite integral_{0}^{+inf} exp(-x) dx  exact=1\n");
    std::printf("    value=%.10g  |err|=%.3e  evals=%zu  status=%s\n",
                semi.value, std::abs(semi.value - exact_semi),
                semi.evaluations, integrator_status_cstr(semi.status));

    // ---- Test 3: integral_{-inf}^{+inf} exp(-x^2) dx = sqrt(pi) --------
    const auto f_gauss = [](double x) { return std::exp(-x*x); };
    AdaptiveResult doubly = integrate_doubly_infinite(f_gauss, opts);
    const double exact_doubly = std::sqrt(M_PI);
    std::printf("\n  Doubly-infinite integral exp(-x^2) dx  exact=sqrt(pi)=%.10g\n", exact_doubly);
    std::printf("    value=%.10g  |err|=%.3e  evals=%zu  status=%s\n\n",
                doubly.value, std::abs(doubly.value - exact_doubly),
                doubly.evaluations, integrator_status_cstr(doubly.status));

    // Write infinite-limit data for plotting.
    {
        FILE* fp = std::fopen("infinite_limits.data", "w");
        if (fp) {
            std::fprintf(fp, "# test exact value abs_err evals\n");
            std::fprintf(fp, "semi_infinite    %.15g %.15g %.3e %zu\n",
                         exact_semi, semi.value, std::abs(semi.value - exact_semi), semi.evaluations);
            std::fprintf(fp, "doubly_infinite  %.15g %.15g %.3e %zu\n",
                         exact_doubly, doubly.value, std::abs(doubly.value - exact_doubly), doubly.evaluations);
            std::fclose(fp);
        }
    }
}

} // namespace pp

// Program entry point.
int main() {
    return pp::run();
}
