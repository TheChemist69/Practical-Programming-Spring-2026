#include "minimization.h"
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace pp;

static double himmelblau(const std::vector<double>& x) {
    double X = x[0];
    double Y = x[1];
    return (X*X + Y - 11.0)*(X*X + Y - 11.0) + (X + Y*Y - 7.0)*(X + Y*Y - 7.0);
}

int main() {
    std::cout << std::fixed << std::setprecision(8);
    
    // Test gradient numerically at (0,0)
    std::vector<double> x = {0.0, 0.0};
    double f0 = himmelblau(x);
    std::cout << "f(0,0) = " << f0 << "\n";
    
    // Try forward difference gradient
    double eps = 1e-8;
    x[0] += eps;
    double fx1 = himmelblau(x);
    x[0] -= eps;
    std::cout << "df/dx (forward) = " << (fx1 - f0) / eps << "\n";
    
    x[1] += eps;
    double fy1 = himmelblau(x);
    x[1] -= eps;
    std::cout << "df/dy (forward) = " << (fy1 - f0) / eps << "\n";
    
    // Try with better starting point
    MinimizerOptions opts;
    opts.acc = 1e-6;
    opts.max_iters = 1000;
    opts.use_central = false;  // use forward differences
    opts.levenberg = 1e-3;      // increase regularization
    
    std::vector<double> x0 = {0.0, 0.0};
    auto res = newton_minimize(himmelblau, x0, opts);
    std::cout << "Result: iterations=" << res.iterations << " converged=" << res.converged
              << " x=[" << res.x[0] << ", " << res.x[1] << "] f=" << res.value << "\n";
    
    return 0;
}
