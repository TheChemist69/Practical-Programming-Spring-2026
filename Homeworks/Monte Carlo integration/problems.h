// "problems.h" header file.
// Canonical benchmark integrals used by this Monte Carlo project.
//
// Each factory returns a fully specified problem:
// - human-readable name,
// - integrand function,
// - integration box [a,b],
// - known exact value (used to compute actual errors).
#pragma once

#include <functional>
#include <string>

#include "vector.h"

namespace pp {

// Self-contained description of one integration test case.
struct MCProblem {
    // Label used in console and data files.
    std::string name;
    // Integrand f(x) to be sampled.
    std::function<double(const vector&)> integrand;
    // Lower bounds of the axis-aligned integration box.
    vector a;
    // Upper bounds of the axis-aligned integration box.
    vector b;
    // Reference value used for accuracy checks.
    double exact = 0.0;
};

// Indicator function of the unit disk on [-1,1]^2, exact integral = pi.
MCProblem unit_circle_area_problem();
// Smooth Gaussian benchmark on [0,1]^2 with closed-form exact value.
MCProblem smooth_gaussian_2d_problem();
// Difficult 3D integrand with endpoint singular behavior near (0,0,0).
MCProblem difficult_singular_3d_problem();

} // namespace pp
