// "spline.h" header file.
// Spline interpolation: linear, quadratic, and cubic splines.
#pragma once

#include "vector.h"

namespace pp {

// -----------------------------------------------------------------
// binsearch -- binary search for the interval containing z.
//
// Given a sorted array x, returns index i such that
//   x[i] <= z <= x[i+1].
// Throws if z is outside [x[0], x[n-1]].
// -----------------------------------------------------------------
int binsearch(const pp::vector& x, NUMBER z);

// -----------------------------------------------------------------
// linterp -- linear spline interpolation.
//
// Given data {x[i], y[i]}, evaluates the piecewise-linear
// interpolant at z.  Uses binsearch to locate the interval.
//
// S_i(z) = y[i] + (y[i+1]-y[i])/(x[i+1]-x[i]) * (z - x[i])
// -----------------------------------------------------------------
NUMBER linterp(const pp::vector& x, const pp::vector& y, NUMBER z);

// -----------------------------------------------------------------
// linterp_integ -- integral of the linear spline from x[0] to z.
//
// Computed analytically as a sum of trapezoids over full intervals
// plus a partial interval:
//
// Int_{x[0]}^{z} S(t) dt = sum_k (y[k]+y[k+1])/2 * h[k]  (full)
//                         + y[i]*u + slope_i*u^2/2          (partial)
// where u = z - x[i], i = binsearch(x, z).
// -----------------------------------------------------------------
NUMBER linterp_integ(const pp::vector& x, const pp::vector& y, NUMBER z);

// -----------------------------------------------------------------
// qspline -- quadratic spline interpolation.
//
// S_i(x) = y[i] + b[i]*(x - x[i]) + c[i]*(x - x[i])^2
//
// Coefficients are computed by averaging forward (c[0]=0) and
// backward (c[m-1]=0) recursions for symmetry:
//
//   c[i]*h[i] + c[i+1]*h[i+1] = p[i+1] - p[i]
//
// where h[i] = x[i+1]-x[i] and p[i] = (y[i+1]-y[i])/h[i].
// -----------------------------------------------------------------
class qspline {
public:
    pp::vector x, y, b, c;

    // Build quadratic spline from data {xs, ys}
    qspline(const pp::vector& xs, const pp::vector& ys);

    // Evaluate spline at z
    NUMBER evaluate(NUMBER z) const;

    // Evaluate first derivative at z
    NUMBER derivative(NUMBER z) const;

    // Evaluate integral from x[0] to z
    NUMBER integral(NUMBER z) const;
};

// -----------------------------------------------------------------
// cspline -- natural cubic spline interpolation.
//
// S_i(x) = y[i] + b[i]*(x-x[i]) + c[i]*(x-x[i])^2 + d[i]*(x-x[i])^3
//
// Natural boundary conditions: S''(x[0]) = S''(x[n-1]) = 0.
// The c[i] coefficients are found by solving a tridiagonal system
// (Thomas algorithm):
//
//   h[i-1]*c[i-1] + 2*(h[i-1]+h[i])*c[i] + h[i]*c[i+1]
//       = 3*(p[i] - p[i-1])
//
// Then:
//   d[i] = (c[i+1] - c[i]) / (3*h[i])
//   b[i] = p[i] - h[i]*(2*c[i] + c[i+1]) / 3
// -----------------------------------------------------------------
class cspline {
public:
    pp::vector x, y, b, c, d;

    // Build natural cubic spline from data {xs, ys}
    cspline(const pp::vector& xs, const pp::vector& ys);

    // Evaluate spline at z
    NUMBER evaluate(NUMBER z) const;

    // Evaluate first derivative at z
    NUMBER derivative(NUMBER z) const;

    // Evaluate integral from x[0] to z
    NUMBER integral(NUMBER z) const;
};

} // namespace pp
