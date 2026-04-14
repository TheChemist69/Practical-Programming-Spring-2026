// "reporting.h" header file.
// Data writer declarations for integration homework outputs.
#pragma once

#include <string>
#include <utility>
#include <vector>

#include "benchmarks.h"
#include "erf_integral.h"

namespace pp {

// Writes one row per benchmark case with pass/fail diagnostics.
void write_benchmark_results(
        const std::string& filename,
        const std::vector<BenchmarkResult>& rows);

// Writes sampled erf curve points with reference and error columns.
void write_erf_curve_data(
        const std::string& filename,
        const std::vector<ErfCurvePoint>& rows);

// Writes tabulated erf points used for visual comparison.
void write_erf_tabulated_data(
        const std::string& filename,
        const std::vector<std::pair<double, double>>& rows);

// Writes erf(1) sweep rows for log-log accuracy plotting.
void write_erf_sweep_data(
        const std::string& filename,
        const std::vector<ErfSweepPoint>& rows);

} // namespace pp
