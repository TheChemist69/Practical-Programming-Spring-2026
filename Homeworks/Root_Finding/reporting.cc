// "reporting.cc" implementation file.
// Writers for Root_Finding homework output tables.

#include "reporting.h"

#include <cmath>
#include <fstream>
#include <iomanip>

namespace pp {

// Writes columns: r, f_numeric, f_exact, |f_numeric - f_exact| for plotting.
void write_wavefunction_data(const std::string& filename, const std::vector<WavefunctionPoint>& rows) {
    std::ofstream out(filename);
    out << std::setprecision(17);
    out << "# r f_numeric f_exact abs_error\n";

    for (const WavefunctionPoint& row : rows) {
        out << row.r << ' '
            << row.f_numeric << ' '
            << row.f_exact << ' '
            << std::abs(row.f_numeric - row.f_exact) << '\n';
    }
}

// Writes a parameter-sweep convergence table with a header line naming
// the varied parameter, then one row per (parameter, energy, error, converged).
void write_convergence_data(const std::string& filename, const std::string& parameter_name, const std::vector<ConvergenceRow>& rows) {
    std::ofstream out(filename);
    out << std::setprecision(17);
    out << "# " << parameter_name << " energy abs_error converged\n";

    for (const ConvergenceRow& row : rows) {
        out << row.parameter << ' '
            << row.energy << ' '
            << row.abs_error << ' '
            << (row.converged ? 1 : 0) << '\n';
    }
}

// Writes one row per (problem, method) pair with all Newton diagnostic counters
// for comparing backtracking vs quadratic-interpolation line search.
void write_line_search_data(const std::string& filename, const std::vector<LineSearchRow>& rows) {
    std::ofstream out(filename);
    out << std::setprecision(17);
    out << "# problem method iterations f_evaluations jacobian_evaluations line_search_iterations residual_norm\n";

    for (const LineSearchRow& row : rows) {
        out << row.problem << ' '
            << row.method << ' '
            << row.iterations << ' '
            << row.f_evaluations << ' '
            << row.jacobian_evaluations << ' '
            << row.line_search_iterations << ' '
            << row.residual_norm << '\n';
    }
}

} // namespace pp
