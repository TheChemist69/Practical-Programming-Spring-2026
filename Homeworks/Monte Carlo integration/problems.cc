// "problems.cc" implementation file.
// Benchmark integrands and exact values for Monte Carlo runs.
//
// These factory functions centralize all problem definitions so the
// experiment driver can reuse the same setup consistently.

#include "problems.h"

#include <cmath>

namespace pp {

MCProblem unit_circle_area_problem() {
    // Integral of indicator(x^2+y^2<=1) over [-1,1]^2 equals pi.
    return MCProblem{
            "unit_circle_area",
            [](const vector& x) {
                const double r2 = x[0] * x[0] + x[1] * x[1];
                return (r2 <= 1.0) ? 1.0 : 0.0;
            },
            vector{-1.0, -1.0},
            vector{1.0, 1.0},
            std::acos(-1.0),
    };
}

MCProblem smooth_gaussian_2d_problem() {
    // Integral of exp(-(x^2+y^2)) on [0,1]^2.
    // In 1D: ∫_0^1 exp(-x^2) dx = sqrt(pi)/2 * erf(1), then square it.
    const double factor = 0.5 * std::sqrt(std::acos(-1.0)) * std::erf(1.0);
    return MCProblem{
            "gaussian_2d",
            [](const vector& x) {
                return std::exp(-(x[0] * x[0] + x[1] * x[1]));
            },
            vector{0.0, 0.0},
            vector{1.0, 1.0},
            factor * factor,
    };
}

MCProblem difficult_singular_3d_problem() {
    // Assignment benchmark:
    // ∫_[0,pi]^3 [1/(pi^3 * (1-cos(x)cos(y)cos(z)))] dxdydz.
    // The denominator approaches zero near the origin, so we protect
    // against exact zero division in floating-point arithmetic.
    const double pi = std::acos(-1.0);
    const double inv_pi3 = 1.0 / (pi * pi * pi);

    return MCProblem{
            "singular_3d_assignment",
            [inv_pi3](const vector& x) {
                const double cprod = std::cos(x[0]) * std::cos(x[1]) * std::cos(x[2]);
                double denom = 1.0 - cprod;
                if (denom < 1e-14) {
                    // Numerical safeguard: keeps the integrand finite.
                    denom = 1e-14;
                }
                return inv_pi3 / denom;
            },
            vector{0.0, 0.0, 0.0},
            vector{pi, pi, pi},
            1.3932039296856768591842462603255,
    };
}

} // namespace pp
