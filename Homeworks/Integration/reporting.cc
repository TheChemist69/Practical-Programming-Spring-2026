// "reporting.cc" implementation file.
// Writers for benchmark and erf data products.

#include "reporting.h"

#include <fstream>
#include <iomanip>

namespace pp {

// Exports benchmark summary table.
void write_benchmark_results(
        const std::string& filename,
        const std::vector<BenchmarkResult>& rows) {
    std::ofstream out(filename);
    out << std::setprecision(17);
    out << "# name a b result exact abs_error est_error evaluations tolerance_goal status pass\n";
    for (const BenchmarkResult& row : rows) {
        out << row.name << ' '
            << row.a << ' '
            << row.b << ' '
            << row.result.value << ' '
            << row.exact << ' '
            << row.abs_error << ' '
            << row.result.estimated_error << ' '
            << row.result.evaluations << ' '
            << row.tolerance_goal << ' '
            << integrator_status_cstr(row.result.status) << ' '
            << (row.passed ? 1 : 0) << '\n';
    }
}

// Exports erf curve samples and pointwise comparison data.
void write_erf_curve_data(
        const std::string& filename,
        const std::vector<ErfCurvePoint>& rows) {
    std::ofstream out(filename);
    out << std::setprecision(17);
    out << "# z erf_numeric erf_reference abs_error est_error evaluations status\n";
    for (const ErfCurvePoint& row : rows) {
        out << row.z << ' '
            << row.value << ' '
            << row.reference << ' '
            << row.abs_error << ' '
            << row.estimated_error << ' '
            << row.evaluations << ' '
            << integrator_status_cstr(row.status) << '\n';
    }
}

// Exports fixed tabulated erf values.
void write_erf_tabulated_data(
        const std::string& filename,
        const std::vector<std::pair<double, double>>& rows) {
    std::ofstream out(filename);
    out << std::setprecision(17);
    out << "# z erf_tabulated\n";
    for (const auto& row : rows) {
        out << row.first << ' ' << row.second << '\n';
    }
}

// Exports erf(1) convergence sweep data.
void write_erf_sweep_data(
        const std::string& filename,
        const std::vector<ErfSweepPoint>& rows) {
    std::ofstream out(filename);
    out << std::setprecision(17);
    out << "# acc erf1_numeric abs_error est_error evaluations status\n";
    for (const ErfSweepPoint& row : rows) {
        out << row.acc << ' '
            << row.value << ' '
            << row.abs_error << ' '
            << row.estimated_error << ' '
            << row.evaluations << ' '
            << integrator_status_cstr(row.status) << '\n';
    }
}

} // namespace pp
