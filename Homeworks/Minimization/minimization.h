// minimization.h
// Newton minimizer with numerical gradient/Hessian and backtracking
#pragma once

#include <vector>
#include <functional>

namespace pp {

struct MinimizerOptions {
    double acc = 1e-6;
    int max_iters = 1000;
    double levenberg = 1e-6;
    bool use_central = true; // central finite differences by default
};

struct MinimizerResult {
    std::vector<double> x;
    int iterations = 0;
    double value = 0.0;
    bool converged = false;
};

// φ : R^n -> double
MinimizerResult newton_minimize(const std::function<double(const std::vector<double>&)>& phi,
                                const std::vector<double>& x0,
                                const MinimizerOptions& opts = {});

}
