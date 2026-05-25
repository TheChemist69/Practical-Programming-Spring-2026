// "benchmarks.h" header file.
// Benchmark integral declarations for adaptive quadrature validation.
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "integrator.h"

namespace pp {

// One benchmark problem: name, integrand, interval, exact result.
struct BenchmarkCase {
    std::string name;
    std::function<double(double)> integrand;
    double a = 0.0;
    double b = 0.0;
    double exact = 0.0;
};

// Numerical outcome for one benchmark case.
struct BenchmarkResult {
    std::string name;
    double a = 0.0;
    double b = 0.0;
    double exact = 0.0;
    AdaptiveResult result;
    double abs_error = 0.0;
    double tolerance_goal = 0.0;
    bool passed = false;
};

// Returns the default assignment benchmark set.
std::vector<BenchmarkCase> default_benchmark_cases();

// Runs all benchmarks with one shared option set.
std::vector<BenchmarkResult> run_benchmarks(const AdaptiveOptions& options);

} // namespace pp
