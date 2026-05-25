// "spline.cc" implementation file.
// Spline interpolation: linear, quadratic, and cubic splines.
// Also includes make_qspline, the functional-programming style variant (Part C).
#include "spline.h"
#include <cmath>
#include <functional>
#include <stdexcept>

namespace pp {

// -----------------------------------------------------------------
// binsearch -- locate interval by bisection.
// Returns i such that x[i] <= z <= x[i+1].
// -----------------------------------------------------------------
int binsearch(const pp::vector& x, NUMBER z) {
    int n = x.size();
    if (z < x[0] || z > x[n - 1]) {
        throw std::runtime_error("binsearch: z out of range");
    }
    int i = 0, j = n - 1;
    while (j - i > 1) {
        int mid = (i + j) / 2;
        if (z > x[mid]) {
            i = mid;
        } else {
            j = mid;
        }
    }
    return i;
}

// -----------------------------------------------------------------
// linterp -- linear interpolation at z.
// S_i(z) = y[i] + dy/dx * (z - x[i])
// -----------------------------------------------------------------
NUMBER linterp(const pp::vector& x, const pp::vector& y, NUMBER z) {
    int i = binsearch(x, z);
    NUMBER dx = x[i + 1] - x[i];
    if (!(dx > 0)) {
        throw std::runtime_error("linterp: dx is not positive");
    }
    NUMBER dy = y[i + 1] - y[i];
    return y[i] + dy / dx * (z - x[i]);
}

// -----------------------------------------------------------------
// linterp_integ -- integral of linear spline from x[0] to z.
// Full intervals: trapezoid (y[k]+y[k+1])/2 * h[k]
// Partial interval: y[i]*u + slope*u^2/2
// -----------------------------------------------------------------
NUMBER linterp_integ(const pp::vector& x, const pp::vector& y, NUMBER z) {
    int i = binsearch(x, z);
    NUMBER sum = 0.0;
    // Sum full intervals [x[k], x[k+1]] for k = 0, ..., i-1
    for (int k = 0; k < i; k++) {
        NUMBER h = x[k + 1] - x[k];
        sum += (y[k] + y[k + 1]) / 2.0 * h;
    }
    // Partial interval [x[i], z]
    NUMBER u = z - x[i];
    NUMBER h = x[i + 1] - x[i];
    NUMBER slope = (y[i + 1] - y[i]) / h;
    // Int_{x[i]}^{z} (y[i] + slope*(t - x[i])) dt = y[i]*u + slope*u^2/2
    sum += y[i] * u + slope * u * u / 2.0;
    return sum;
}

// -----------------------------------------------------------------
// qspline -- quadratic spline constructor.
//
// Forward recursion (c_f[0] = 0):
//   c_f[i+1] = (p[i+1] - p[i] - c_f[i]*h[i]) / h[i+1]
//
// Backward recursion (c_b[m-1] = 0):
//   c_b[i]   = (p[i+1] - p[i] - c_b[i+1]*h[i+1]) / h[i]
//
// Average: c[i] = (c_f[i] + c_b[i]) / 2
// Then:    b[i] = p[i] - c[i]*h[i]
// -----------------------------------------------------------------
qspline::qspline(const pp::vector& xs, const pp::vector& ys) : x(xs), y(ys) {
    int n = x.size();
    int m = n - 1;  // number of intervals
    if (n < 2) {
        throw std::runtime_error("qspline: need at least 2 data points");
    }

    // Interval widths h[i] and slopes p[i]
    pp::vector h(m), p(m);
    for (int i = 0; i < m; i++) {
        h[i] = x[i + 1] - x[i];
        if (!(h[i] > 0)) {
            throw std::runtime_error("qspline: x must be strictly increasing");
        }
        p[i] = (y[i + 1] - y[i]) / h[i];
    }

    // Forward recursion: c_f[0] = 0
    pp::vector cf(m);
    cf[0] = 0.0;
    for (int i = 0; i < m - 1; i++) {
        cf[i + 1] = (p[i + 1] - p[i] - cf[i] * h[i]) / h[i + 1];
    }

    // Backward recursion: c_b[m-1] = 0
    pp::vector cb(m);
    cb[m - 1] = 0.0;
    for (int i = m - 2; i >= 0; i--) {
        cb[i] = (p[i + 1] - p[i] - cb[i + 1] * h[i + 1]) / h[i];
    }

    // Average and compute b
    c.resize(m);
    b.resize(m);
    for (int i = 0; i < m; i++) {
        c[i] = (cf[i] + cb[i]) / 2.0;
        b[i] = p[i] - c[i] * h[i];
    }
}

// S_i(z) = y[i] + b[i]*u + c[i]*u^2,  u = z - x[i]
NUMBER qspline::evaluate(NUMBER z) const {
    int i = binsearch(x, z);
    NUMBER u = z - x[i];
    return y[i] + u * (b[i] + c[i] * u);
}

// S'_i(z) = b[i] + 2*c[i]*u
NUMBER qspline::derivative(NUMBER z) const {
    int i = binsearch(x, z);
    NUMBER u = z - x[i];
    return b[i] + 2.0 * c[i] * u;
}

// Int_{x[0]}^{z} S(t) dt
// Full interval: y[k]*h + b[k]*h^2/2 + c[k]*h^3/3
// Partial:       y[i]*u + b[i]*u^2/2 + c[i]*u^3/3
NUMBER qspline::integral(NUMBER z) const {
    int i = binsearch(x, z);
    NUMBER sum = 0.0;
    for (int k = 0; k < i; k++) {
        NUMBER h = x[k + 1] - x[k];
        sum += y[k] * h + b[k] * h * h / 2.0 + c[k] * h * h * h / 3.0;
    }
    NUMBER u = z - x[i];
    sum += y[i] * u + b[i] * u * u / 2.0 + c[i] * u * u * u / 3.0;
    return sum;
}

// -----------------------------------------------------------------
// cspline -- natural cubic spline constructor.
//
// Solves the tridiagonal system for c[1], ..., c[n-2]:
//   h[i-1]*c[i-1] + 2*(h[i-1]+h[i])*c[i] + h[i]*c[i+1]
//       = 3*(p[i] - p[i-1])
// with c[0] = c[n-1] = 0  (natural boundary conditions).
//
// Then:
//   d[i] = (c[i+1] - c[i]) / (3*h[i])
//   b[i] = p[i] - h[i]*(2*c[i] + c[i+1]) / 3
// -----------------------------------------------------------------
cspline::cspline(const pp::vector& xs, const pp::vector& ys) : x(xs), y(ys) {
    int n = x.size();
    int m = n - 1;  // number of intervals
    if (n < 2) {
        throw std::runtime_error("cspline: need at least 2 data points");
    }

    // Interval widths and slopes
    pp::vector h(m), p(m);
    for (int i = 0; i < m; i++) {
        h[i] = x[i + 1] - x[i];
        if (!(h[i] > 0)) {
            throw std::runtime_error("cspline: x must be strictly increasing");
        }
        p[i] = (y[i + 1] - y[i]) / h[i];
    }

    // Natural boundary: c[0] = 0, c[n-1] = 0
    c.resize(n);
    c[0] = 0.0;
    c[n - 1] = 0.0;

    int sz = n - 2;  // number of interior unknowns c[1]...c[n-2]
    if (sz > 0) {
        // Build tridiagonal system
        // Sub-diagonal a[], main diagonal diag[], super-diagonal e[], rhs[]
        pp::vector a(sz), diag(sz), e(sz), rhs(sz);
        for (int j = 0; j < sz; j++) {
            int idx = j + 1;  // c-index: c[idx]
            diag[j] = 2.0 * (h[idx - 1] + h[idx]);
            rhs[j]  = 3.0 * (p[idx] - p[idx - 1]);
            if (j > 0) {
                a[j] = h[idx - 1];
            }
            if (j < sz - 1) {
                e[j] = h[idx];
            }
        }

        // Thomas algorithm -- forward elimination
        for (int j = 1; j < sz; j++) {
            NUMBER w = a[j] / diag[j - 1];
            diag[j] -= w * e[j - 1];
            rhs[j]  -= w * rhs[j - 1];
        }

        // Thomas algorithm -- back substitution
        pp::vector csol(sz);
        csol[sz - 1] = rhs[sz - 1] / diag[sz - 1];
        for (int j = sz - 2; j >= 0; j--) {
            csol[j] = (rhs[j] - e[j] * csol[j + 1]) / diag[j];
        }

        // Copy into full c array
        for (int j = 0; j < sz; j++) {
            c[j + 1] = csol[j];
        }
    }

    // Compute b[i] and d[i] for each interval
    b.resize(m);
    d.resize(m);
    for (int i = 0; i < m; i++) {
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        b[i] = p[i] - h[i] * (2.0 * c[i] + c[i + 1]) / 3.0;
    }
}

// S_i(z) = y[i] + b[i]*u + c[i]*u^2 + d[i]*u^3
NUMBER cspline::evaluate(NUMBER z) const {
    int i = binsearch(x, z);
    NUMBER u = z - x[i];
    return y[i] + u * (b[i] + u * (c[i] + u * d[i]));
}

// S'_i(z) = b[i] + 2*c[i]*u + 3*d[i]*u^2
NUMBER cspline::derivative(NUMBER z) const {
    int i = binsearch(x, z);
    NUMBER u = z - x[i];
    return b[i] + u * (2.0 * c[i] + 3.0 * d[i] * u);
}

// Int_{x[0]}^{z} S(t) dt
// Full interval: y[k]*h + b[k]*h^2/2 + c[k]*h^3/3 + d[k]*h^4/4
// Partial:       y[i]*u + b[i]*u^2/2 + c[i]*u^3/3 + d[i]*u^4/4
NUMBER cspline::integral(NUMBER z) const {
    int i = binsearch(x, z);
    NUMBER sum = 0.0;
    for (int k = 0; k < i; k++) {
        NUMBER hk = x[k + 1] - x[k];
        sum += y[k] * hk
             + b[k] * hk * hk / 2.0
             + c[k] * hk * hk * hk / 3.0
             + d[k] * hk * hk * hk * hk / 4.0;
    }
    NUMBER u = z - x[i];
    sum += y[i] * u
         + b[i] * u * u / 2.0
         + c[i] * u * u * u / 3.0
         + d[i] * u * u * u * u / 4.0;
    return sum;
}

// -----------------------------------------------------------------
// make_qspline -- functional-programming style (Part C).
//
// Computes b and c exactly as the qspline constructor does, then
// move-captures x, y, b, c into a lambda by value.  The returned
// std::function<double(double)> is the spline evaluator.
// Returning a closure rather than an object is the functional style:
// the caller gets a plain callable with no visible internal state.
// -----------------------------------------------------------------
std::function<double(double)> make_qspline(const pp::vector& xs, const pp::vector& ys) {
    int n = xs.size();
    int m = n - 1;
    if (n < 2) throw std::runtime_error("make_qspline: need at least 2 data points");

    pp::vector h(m), p(m);
    for (int i = 0; i < m; i++) {
        h[i] = xs[i+1] - xs[i];
        if (!(h[i] > 0)) throw std::runtime_error("make_qspline: x must be strictly increasing");
        p[i] = (ys[i+1] - ys[i]) / h[i];
    }

    // Forward recursion (c_f[0] = 0)
    pp::vector cf(m);
    cf[0] = 0.0;
    for (int i = 0; i < m-1; i++)
        cf[i+1] = (p[i+1] - p[i] - cf[i]*h[i]) / h[i+1];

    // Backward recursion (c_b[m-1] = 0)
    pp::vector cb(m);
    cb[m-1] = 0.0;
    for (int i = m-2; i >= 0; i--)
        cb[i] = (p[i+1] - p[i] - cb[i+1]*h[i+1]) / h[i];

    // Averaged coefficients
    pp::vector b(m), c(m);
    for (int i = 0; i < m; i++) {
        c[i] = (cf[i] + cb[i]) / 2.0;
        b[i] = p[i] - c[i]*h[i];
    }

    // Capture x, y, b, c by value (move-copy) into the closure.
    // The closure is self-contained: no reference to local variables.
    pp::vector x_cap = xs;
    pp::vector y_cap = ys;
    pp::vector b_cap = std::move(b);
    pp::vector c_cap = std::move(c);

    return [x_cap = std::move(x_cap),
            y_cap = std::move(y_cap),
            b_cap = std::move(b_cap),
            c_cap = std::move(c_cap)](double z) -> double {
        int idx = binsearch(x_cap, z);
        double u = z - x_cap[idx];
        return y_cap[idx] + u * (b_cap[idx] + c_cap[idx] * u);
    };
}

} // namespace pp
