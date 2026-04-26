// "reporting.h" header file.
// Data-file writers for plots and diagnostics.
#pragma once

#include <string>
#include <vector>

#include "hydrogen.h"

namespace pp {

struct LineSearchRow {
    std::string problem;
    std::string method;
    int iterations = 0;
    int f_evaluations = 0;
    int jacobian_evaluations = 0;
    int line_search_iterations = 0;
    double residual_norm = 0.0;
};

void write_wavefunction_data(const std::string& filename, const std::vector<WavefunctionPoint>& rows);
void write_convergence_data(const std::string& filename, const std::string& parameter_name, const std::vector<ConvergenceRow>& rows);
void write_line_search_data(const std::string& filename, const std::vector<LineSearchRow>& rows);

} // namespace pp
