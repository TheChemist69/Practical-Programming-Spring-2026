// "reporting.h" header file.
// Data schemas and file writers for all experiment outputs.
//
// Each struct mirrors one .data file row layout so plotting scripts can
// consume files without additional transformations.
#pragma once

#include <string>
#include <vector>

namespace pp {

// One row for the quick two-dimensional sanity examples.
struct TwoDExampleRow {
    // Algorithm label ("plain", etc.).
    std::string method;
    // Problem label (matches MCProblem::name).
    std::string problem;
    // Number of samples used.
    int n = 0;
    // Estimated integral.
    double value = 0.0;
    // Exact reference value.
    double exact = 0.0;
    // |value - exact|.
    double actual_error = 0.0;
    // Method-reported error estimate.
    double estimated_error = 0.0;
};

// One row for plain Monte Carlo convergence study.
struct PlainScalingRow {
    int n = 0;
    double estimated_error = 0.0;
    double actual_error = 0.0;
    double value = 0.0;
    double exact = 0.0;
    // Error * sqrt(N) helps inspect asymptotic constancy.
    double estimated_times_sqrt_n = 0.0;
    double actual_times_sqrt_n = 0.0;
};

// One row comparing plain pseudo-random vs quasi-random MC.
struct QuasiScalingRow {
    int n = 0;
    double actual_error_plain = 0.0;
    double actual_error_quasi = 0.0;
    double estimated_error_plain = 0.0;
    double estimated_error_quasi = 0.0;
    double value_plain = 0.0;
    double value_quasi = 0.0;
    double exact = 0.0;
};

// One row comparing plain MC and stratified MC.
struct StratifiedScalingRow {
    int n = 0;
    double actual_error_plain = 0.0;
    double actual_error_stratified = 0.0;
    double estimated_error_plain = 0.0;
    double estimated_error_stratified = 0.0;
    // Number of integrand evaluations effectively spent by each method.
    double calls_plain = 0.0;
    double calls_stratified = 0.0;
    double exact = 0.0;
};

// One row for the difficult singular benchmark where all methods are shown.
struct DifficultIntegralRow {
    int n = 0;
    double value_plain = 0.0;
    double value_quasi = 0.0;
    double value_stratified = 0.0;
    double actual_error_plain = 0.0;
    double actual_error_quasi = 0.0;
    double actual_error_stratified = 0.0;
    double exact = 0.0;
};

// Write the corresponding table to disk with a fixed column order.
void write_two_dim_examples(const std::string& filename, const std::vector<TwoDExampleRow>& rows);
void write_plain_scaling(const std::string& filename, const std::vector<PlainScalingRow>& rows);
void write_quasi_scaling(const std::string& filename, const std::vector<QuasiScalingRow>& rows);
void write_stratified_scaling(const std::string& filename, const std::vector<StratifiedScalingRow>& rows);
void write_difficult_integral(const std::string& filename, const std::vector<DifficultIntegralRow>& rows);

} // namespace pp
