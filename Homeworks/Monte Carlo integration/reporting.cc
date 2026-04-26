// "reporting.cc" implementation file.
// Writers for Monte Carlo experiment data files.
//
// Each writer emits:
// 1) one header line documenting column semantics,
// 2) one line per row in the same order expected by gnuplot scripts.

#include "reporting.h"

#include <fstream>
#include <iomanip>

namespace pp {

void write_two_dim_examples(const std::string& filename, const std::vector<TwoDExampleRow>& rows) {
    // Quick two-problem summary table used for sanity checking.
    std::ofstream out(filename);
    out << std::setprecision(17);
    out << "# method problem N value exact actual_error estimated_error\n";
    for (const TwoDExampleRow& row : rows) {
        out << row.method << ' '
            << row.problem << ' '
            << row.n << ' '
            << row.value << ' '
            << row.exact << ' '
            << row.actual_error << ' '
            << row.estimated_error << '\n';
    }
}

void write_plain_scaling(const std::string& filename, const std::vector<PlainScalingRow>& rows) {
    // Convergence table for plain pseudo-random Monte Carlo.
    std::ofstream out(filename);
    out << std::setprecision(17);
    out << "# N estimated_error actual_error value exact estimated_error_sqrtN actual_error_sqrtN\n";
    for (const PlainScalingRow& row : rows) {
        out << row.n << ' '
            << row.estimated_error << ' '
            << row.actual_error << ' '
            << row.value << ' '
            << row.exact << ' '
            << row.estimated_times_sqrt_n << ' '
            << row.actual_times_sqrt_n << '\n';
    }
}

void write_quasi_scaling(const std::string& filename, const std::vector<QuasiScalingRow>& rows) {
    // Side-by-side table comparing quasi-random and plain MC.
    std::ofstream out(filename);
    out << std::setprecision(17);
    out << "# N actual_plain actual_quasi estimated_plain estimated_quasi value_plain value_quasi exact\n";
    for (const QuasiScalingRow& row : rows) {
        out << row.n << ' '
            << row.actual_error_plain << ' '
            << row.actual_error_quasi << ' '
            << row.estimated_error_plain << ' '
            << row.estimated_error_quasi << ' '
            << row.value_plain << ' '
            << row.value_quasi << ' '
            << row.exact << '\n';
    }
}

void write_stratified_scaling(const std::string& filename, const std::vector<StratifiedScalingRow>& rows) {
    // Comparison of stratified and plain MC including evaluation counts.
    std::ofstream out(filename);
    out << std::setprecision(17);
    out << "# N actual_plain actual_stratified estimated_plain estimated_stratified calls_plain calls_stratified exact\n";
    for (const StratifiedScalingRow& row : rows) {
        out << row.n << ' '
            << row.actual_error_plain << ' '
            << row.actual_error_stratified << ' '
            << row.estimated_error_plain << ' '
            << row.estimated_error_stratified << ' '
            << row.calls_plain << ' '
            << row.calls_stratified << ' '
            << row.exact << '\n';
    }
}

void write_difficult_integral(const std::string& filename, const std::vector<DifficultIntegralRow>& rows) {
    // Method comparison on the difficult singular benchmark.
    std::ofstream out(filename);
    out << std::setprecision(17);
    out << "# N value_plain value_quasi value_stratified actual_plain actual_quasi actual_stratified exact\n";
    for (const DifficultIntegralRow& row : rows) {
        out << row.n << ' '
            << row.value_plain << ' '
            << row.value_quasi << ' '
            << row.value_stratified << ' '
            << row.actual_error_plain << ' '
            << row.actual_error_quasi << ' '
            << row.actual_error_stratified << ' '
            << row.exact << '\n';
    }
}

} // namespace pp
