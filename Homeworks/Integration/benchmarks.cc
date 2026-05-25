// "benchmarks.cc" implementation file.
// Benchmark integral definitions and runner.

#include "benchmarks.h"

#include <cmath>

namespace pp {

// Core assignment benchmarks, including endpoint-singular examples.
std::vector<BenchmarkCase> default_benchmark_cases() {
    const double pi = std::acos(-1.0);
    return {
        {
            "int_0_1_sqrt_x",
            [](double x) { return std::sqrt(x); },
            0.0,
            1.0,
            2.0 / 3.0,
        },
        {
            "int_0_1_inv_sqrt_x",
            [](double x) { return 1.0 / std::sqrt(x); },
            0.0,
            1.0,
            2.0,
        },
        {
            "int_0_1_sqrt_1_minus_x2",
            [](double x) { return std::sqrt(1.0 - x * x); },
            0.0,
            1.0,
            pi / 4.0,
        },
        {
            "int_0_1_log_x_over_sqrt_x",
            [](double x) { return std::log(x) / std::sqrt(x); },
            0.0,
            1.0,
            -4.0,
        },
    };
}

// Evaluates each case and records absolute error and pass/fail.
std::vector<BenchmarkResult> run_benchmarks(const AdaptiveOptions& options) {
    const std::vector<BenchmarkCase> cases = default_benchmark_cases();
    std::vector<BenchmarkResult> rows;
    rows.reserve(cases.size());

    for (const BenchmarkCase& c : cases) {
        BenchmarkResult row;
        row.name = c.name;
        row.a = c.a;
        row.b = c.b;
        row.exact = c.exact;
        row.result = integrate_open4_adaptive(c.integrand, c.a, c.b, options);
        row.abs_error = std::abs(row.result.value - c.exact);
        row.tolerance_goal = options.acc + options.eps * std::abs(row.result.value);
        row.passed = (row.result.status == IntegratorStatus::success) && (row.abs_error <= row.tolerance_goal);
        rows.push_back(row);
    }

    return rows;
}

} // namespace pp
