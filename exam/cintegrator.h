// "cintegrator.h" header file.
// Recursive open-4-point adaptive integrator (book Table 7.3) for COMPLEX-valued
// integrands, plus straight-line contour helpers.  Everything works in
// std::complex<double>: a real integral is one whose imaginary part is zero,
// and the caller takes the real part.  (See cintegrator.cc for the rule.)
#pragma once

#include <complex>
#include <cstddef>
#include <functional>
#include <vector>

namespace pp {

using cplx = std::complex<double>;

// Result returned by every integrator entry point.
struct CResult {
    cplx value{};                 // estimated integral
    std::size_t evaluations = 0;  // total integrand calls
};

// Tunable accuracy goals and recursion guard.
struct COptions {
    double acc = 1e-8;     // absolute accuracy goal (delta)
    double eps = 1e-8;     // relative accuracy goal (epsilon)
    int max_depth = 1000;  // recursion guard (smooth integrands need only tens of levels)
};

// Adaptive integral of a complex-valued function g over the real interval [a,b].
// This is the engine; the contour helpers below build the appropriate g.
CResult adaptive_unit(const std::function<cplx(double)>& g,
                      double a, double b, const COptions& opt = {});

// Complex line integral of f along the straight segment from a to b:
//   integral = (b - a) * \int_0^1 f(a + t(b-a)) dt .
CResult integrate_segment(const std::function<cplx(cplx)>& f,
                          cplx a, cplx b, const COptions& opt = {});

// Complex integral around a closed polygon through the given vertices, in order,
// with the closing edge (last -> first) added automatically.
CResult integrate_polygon(const std::function<cplx(cplx)>& f,
                          const std::vector<cplx>& vertices,
                          const COptions& opt = {});

}  // namespace pp
