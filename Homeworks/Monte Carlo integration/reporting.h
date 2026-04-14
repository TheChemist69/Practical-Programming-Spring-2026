// "reporting.h" header file.
// Data row declarations and writer interfaces for Monte Carlo experiments.
#pragma once

#include <string>
#include <vector>

namespace pp {

struct TwoDExampleRow {
    std::string method;
    std::string problem;
    int n = 0;
    double value = 0.0;
    double exact = 0.0;
    double actual_error = 0.0;
    double estimated_error = 0.0;
};

struct PlainScalingRow {
    int n = 0;
    double estimated_error = 0.0;
    double actual_error = 0.0;
    double value = 0.0;
    double exact = 0.0;
    double estimated_times_sqrt_n = 0.0;
    double actual_times_sqrt_n = 0.0;
};

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

struct StratifiedScalingRow {
    int n = 0;
    double actual_error_plain = 0.0;
    double actual_error_stratified = 0.0;
    double estimated_error_plain = 0.0;
    double estimated_error_stratified = 0.0;
    double calls_plain = 0.0;
    double calls_stratified = 0.0;
    double exact = 0.0;
};

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

void write_two_dim_examples(const std::string& filename, const std::vector<TwoDExampleRow>& rows);
void write_plain_scaling(const std::string& filename, const std::vector<PlainScalingRow>& rows);
void write_quasi_scaling(const std::string& filename, const std::vector<QuasiScalingRow>& rows);
void write_stratified_scaling(const std::string& filename, const std::vector<StratifiedScalingRow>& rows);
void write_difficult_integral(const std::string& filename, const std::vector<DifficultIntegralRow>& rows);

} // namespace pp
