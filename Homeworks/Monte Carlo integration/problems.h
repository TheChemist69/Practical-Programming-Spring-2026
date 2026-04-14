// "problems.h" header file.
// Benchmark integrands and exact values for Monte Carlo homework.
#pragma once

#include <functional>
#include <string>

#include "vector.h"

namespace pp {

struct MCProblem {
    std::string name;
    std::function<double(const vector&)> integrand;
    vector a;
    vector b;
    double exact = 0.0;
};

MCProblem unit_circle_area_problem();
MCProblem smooth_gaussian_2d_problem();
MCProblem difficult_singular_3d_problem();

} // namespace pp
