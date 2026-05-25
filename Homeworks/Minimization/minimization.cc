// minimization.cc
// Newton minimizer with numerical gradient/Hessian and backtracking line search.
//
// The minimizer builds a local quadratic model at each step:
//   phi(x + dx) ≈ phi(x) + g^T dx + 0.5 dx^T H dx
// and solves (H + lambda*I) dx = -g for the Newton step.
// A backtracking line search then accepts the step if phi decreases.

#include "minimization.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>

namespace pp {

using vec = std::vector<double>;
using mat = std::vector<std::vector<double>>;

// Euclidean norm of a vector.
static double norm2(const vec& v) {
    double s = 0.0;
    for (double a : v) s += a*a;
    return std::sqrt(s);
}

// Solves Ax = b by Gaussian elimination with partial pivoting.
// Returns false if the system is (near-)singular.
static bool solve_linear_system(mat A, vec b, vec& x) {
    int n = (int)b.size();
    x.assign(n, 0.0);
    const double EPS = 1e-16;
    for (int i = 0; i < n; ++i) {
        int piv = i;
        for (int k = i+1; k < n; ++k) if (std::abs(A[k][i]) > std::abs(A[piv][i])) piv = k;
        if (std::abs(A[piv][i]) < EPS) return false;
        if (piv != i) { std::swap(A[piv], A[i]); std::swap(b[piv], b[i]); }
        double diag = A[i][i];
        for (int j = i; j < n; ++j) A[i][j] /= diag;
        b[i] /= diag;
        for (int r = 0; r < n; ++r) if (r != i) {
            double f = A[r][i];
            for (int c = i; c < n; ++c) A[r][c] -= f * A[i][c];
            b[r] -= f * b[i];
        }
    }
    for (int i = 0; i < n; ++i) x[i] = b[i];
    return true;
}

// Computes the gradient of phi at x using finite differences.
// Forward differences: O(h) accuracy, optimal h ~ eps^(1/2) ~ 2^(-26).
// Central differences: O(h^2) accuracy, optimal h ~ eps^(1/3) ~ 2^(-17).
// Using the correct step sizes is essential for Part C's comparison to be meaningful.
static vec numerical_gradient(const std::function<double(const vec&)>& phi, const vec& x, bool central) {
    int n = (int)x.size();
    vec g(n);
    double fx = phi(x);
    vec xt = x;
    for (int i = 0; i < n; ++i) {
        // Forward: optimal h = eps^(1/2) * (1 + |x_i|),  eps = 2^(-52)  =>  h ~ 2^(-26)
        // Central: optimal h = eps^(1/3) * (1 + |x_i|),  eps = 2^(-52)  =>  h ~ 2^(-17)
        double dx = (1.0 + std::abs(x[i])) * std::pow(2.0, central ? -17.0 : -26.0);
        if (central) {
            xt[i] = x[i] + dx;
            double f1 = phi(xt);
            xt[i] = x[i] - dx;
            double f2 = phi(xt);
            g[i] = (f1 - f2) / (2.0*dx);
            xt[i] = x[i];
        } else {
            xt[i] = x[i] + dx;
            double f1 = phi(xt);
            g[i] = (f1 - fx) / dx;
            xt[i] = x[i];
        }
    }
    return g;
}

// Computes the Hessian of phi at x by differencing the gradient.
// Uses the same forward/central flag for consistency with the gradient.
static mat numerical_hessian(const std::function<double(const vec&)>& phi, const vec& x, bool central) {
    int n = (int)x.size();
    mat H(n, vec(n, 0.0));
    vec g0 = numerical_gradient(phi, x, central);
    vec xt = x;
    for (int j = 0; j < n; ++j) {
        // Step for Hessian differencing: cbrt(eps) is optimal for central, 2^(-13) for forward.
        double dxj;
        if (central) dxj = std::cbrt(std::numeric_limits<double>::epsilon()) * (1.0+std::abs(x[j]));
        else         dxj = (1.0 + std::abs(x[j])) * std::pow(2.0, -13.0);
        xt[j] = x[j] + dxj;
        vec g1 = numerical_gradient(phi, xt, central);
        for (int i = 0; i < n; ++i) H[i][j] = (g1[i] - g0[i]) / dxj;
        xt[j] = x[j];
    }
    return H;
}

// Newton minimizer with Levenberg regularization and backtracking line search.
// At each iteration, solves (H + levenberg*I) dx = -g, then backtracks until phi decreases.
// The levenberg parameter should be small (1e-6) so the method stays close to Newton's method.
// Larger values damp towards gradient descent and distort iteration-count comparisons.
MinimizerResult newton_minimize(const std::function<double(const vec&)>& phi,
                                const vec& x0,
                                const MinimizerOptions& opts) {
    vec x = x0;
    MinimizerResult res;
    res.x = x;
    for (int iter = 0; iter < opts.max_iters; ++iter) {
        vec g = numerical_gradient(phi, x, opts.use_central);
        double gnorm = norm2(g);
        res.iterations = iter;
        res.value = phi(x);
        if (gnorm < opts.acc) { res.converged = true; res.x = x; return res; }
        mat H = numerical_hessian(phi, x, opts.use_central);
        int n = (int)x.size();
        // Small Levenberg shift keeps H positive-definite near saddle points.
        for (int i = 0; i < n; ++i) H[i][i] += opts.levenberg;
        vec b(n);
        for (int i = 0; i < n; ++i) b[i] = -g[i];
        vec dx;
        bool solved = solve_linear_system(H, b, dx);
        if (!solved) {
            // Singular H: add a larger shift and retry once.
            for (int i = 0; i < n; ++i) H[i][i] += 1e-3;
            if (!solve_linear_system(H, b, dx)) { res.converged = false; res.x = x; return res; }
        }
        // Backtracking: halve the step until phi decreases.
        double lambda = 1.0;
        double fx = phi(x);
        vec xnew(n);
        while (lambda >= 1.0/1024.0) {
            for (int i = 0; i < n; ++i) xnew[i] = x[i] + lambda * dx[i];
            double fnew = phi(xnew);
            if (fnew < fx) { x = xnew; break; }
            lambda *= 0.5;
        }
        res.x = x;
    }
    res.converged = false;
    res.value = phi(x);
    return res;
}

}
