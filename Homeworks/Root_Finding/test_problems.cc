// "test_problems.cc" implementation file.
// Benchmark equations and analytic gradients for root finding.

#include "test_problems.h"

#include <cassert>

namespace pp {

VectorFunction simple_scalar_equation() {
    // Root at x = sqrt(2).
    return [](const vector& x) {
        assert(x.size() == 1);
        return vector{x[0] * x[0] - 2.0};
    };
}

VectorFunction simple_two_dimensional_system() {
    // Root at (x,y) = (1,2).
    return [](const vector& x) {
        assert(x.size() == 2);
        const double xx = x[0];
        const double yy = x[1];
        return vector{xx + yy - 3.0, xx * xx + yy * yy - 5.0};
    };
}

double rosenbrock_value(const vector& x) {
    assert(x.size() == 2);
    const double xx = x[0];
    const double yy = x[1];
    const double a = 1.0 - xx;
    const double b = yy - xx * xx;
    return a * a + 100.0 * b * b;
}

VectorFunction rosenbrock_gradient() {
    return [](const vector& x) {
        assert(x.size() == 2);
        const double xx = x[0];
        const double yy = x[1];
        const double gx = -2.0 * (1.0 - xx) - 400.0 * xx * (yy - xx * xx);
        const double gy = 200.0 * (yy - xx * xx);
        return vector{gx, gy};
    };
}

double himmelblau_value(const vector& x) {
    assert(x.size() == 2);
    const double xx = x[0];
    const double yy = x[1];
    const double u = xx * xx + yy - 11.0;
    const double v = xx + yy * yy - 7.0;
    return u * u + v * v;
}

VectorFunction himmelblau_gradient() {
    return [](const vector& x) {
        assert(x.size() == 2);
        const double xx = x[0];
        const double yy = x[1];
        const double u = xx * xx + yy - 11.0;
        const double v = xx + yy * yy - 7.0;
        const double gx = 4.0 * xx * u + 2.0 * v;
        const double gy = 2.0 * u + 4.0 * yy * v;
        return vector{gx, gy};
    };
}

} // namespace pp
