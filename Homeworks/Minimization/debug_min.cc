#include <vector>
#include <functional>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>

using vec = std::vector<double>;
using mat = std::vector<std::vector<double>>;

static double norm2(const vec& v) {
    double s = 0.0;
    for (double a : v) s += a*a;
    return std::sqrt(s);
}

static vec numerical_gradient_simple(const std::function<double(const vec&)>& phi, const vec& x) {
    int n = (int)x.size();
    vec g(n);
    double fx = phi(x);
    vec xt = x;
    double eps = std::numeric_limits<double>::epsilon();
    double delta = std::sqrt(eps);
    for (int i = 0; i < n; ++i) {
        double dx = delta * (1.0 + std::abs(x[i]));
        xt[i] = x[i] + dx;
        double f1 = phi(xt);
        g[i] = (f1 - fx) / dx;
        xt[i] = x[i];
    }
    return g;
}

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

static double himmelblau(const std::vector<double>& x) {
    double X = x[0];
    double Y = x[1];
    return (X*X + Y - 11.0)*(X*X + Y - 11.0) + (X + Y*Y - 7.0)*(X + Y*Y - 7.0);
}

int main() {
    std::cout << std::fixed << std::setprecision(8);
    
    vec x = {0.0, 0.0};
    
    // Test one iteration
    std::cout << "Starting at x=[" << x[0] << ", " << x[1] << "]\n";
    std::cout << "f(x) = " << himmelblau(x) << "\n";
    
    vec g = numerical_gradient_simple(himmelblau, x);
    std::cout << "gradient = [" << g[0] << ", " << g[1] << "]\n";
    std::cout << "grad norm = " << norm2(g) << "\n";
    
    // Approximate Hessian with finite differences
    mat H(2, vec(2, 0.0));
    vec xt = x;
    double eps = std::numeric_limits<double>::epsilon();
    double delta = std::pow(eps, 0.25);
    
    for (int j = 0; j < 2; ++j) {
        double dxj = delta * (1.0 + std::abs(x[j]));
        xt[j] = x[j] + dxj;
        vec g1 = numerical_gradient_simple(himmelblau, xt);
        for (int i = 0; i < 2; ++i) {
            H[i][j] = (g1[i] - g[i]) / dxj;
        }
        xt[j] = x[j];
    }
    
    std::cout << "Hessian = [[" << H[0][0] << ", " << H[0][1] << "]\n";
    std::cout << "           [" << H[1][0] << ", " << H[1][1] << "]]\n";
    
    // Add Levenberg
    H[0][0] += 1e-3;
    H[1][1] += 1e-3;
    std::cout << "After Levenberg: [[" << H[0][0] << ", " << H[0][1] << "]\n";
    std::cout << "                  [" << H[1][0] << ", " << H[1][1] << "]]\n";
    
    // Solve H*dx = -g
    vec b(2);
    b[0] = -g[0];
    b[1] = -g[1];
    vec dx;
    bool ok = solve_linear_system(H, b, dx);
    std::cout << "Linear solve ok: " << ok << "\n";
    std::cout << "dx = [" << dx[0] << ", " << dx[1] << "]\n";
    
    // Test linesearch
    double fx = himmelblau(x);
    vec xnew(2);
    double lambda = 1.0;
    while (lambda >= 1.0/1024.0) {
        xnew[0] = x[0] + lambda * dx[0];
        xnew[1] = x[1] + lambda * dx[1];
        double fnew = himmelblau(xnew);
        std::cout << "lambda=" << lambda << " fnew=" << fnew << " < " << fx << " ? " << (fnew < fx ? "yes" : "no") << "\n";
        if (fnew < fx) break;
        lambda *= 0.5;
    }
    
    return 0;
}
