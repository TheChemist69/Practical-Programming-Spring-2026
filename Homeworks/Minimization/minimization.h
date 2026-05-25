// "minimization.h" header file.
// Newton minimizer with numerical gradient/Hessian and backtracking line search.
//
// Algorithm:
//   At each iteration the gradient g = nabla*phi(x) and Hessian H = nabla^2*phi(x)
//   are approximated by finite differences.  The Newton step dx = -(H + lm*I)^{-1} g
//   uses Levenberg regularization (lm*I) to keep H positive-definite near flat regions.
//   A backtracking line search with the Armijo condition ensures descent.
//   Convergence is declared when ||g|| < acc.
#pragma once

#include <vector>
#include <functional>

namespace pp {

// Tunable parameters for the Newton minimizer.
struct MinimizerOptions {
    double acc = 1e-6;       // convergence threshold: stop when ||grad phi(x)|| < acc
    int max_iters = 1000;    // maximum number of Newton iterations
    double levenberg = 1e-6; // Levenberg regularization added to Hessian diagonal (H + lm*I);
                             // keeps the system positive-definite without distorting Newton character
    bool use_central = true; // use central differences for gradient/Hessian (more accurate than forward)
};

// Output of a newton_minimize call.
struct MinimizerResult {
    std::vector<double> x;  // best point found (minimizer or last iterate)
    int iterations = 0;     // number of Newton steps taken
    double value = 0.0;     // phi(x) at the returned point
    bool converged = false; // true if ||grad phi(x)|| < acc was reached
};

// Finds a local minimum of phi : R^n -> R starting from x0.
// Uses numerical gradient/Hessian via finite differences and backtracking line search.
MinimizerResult newton_minimize(const std::function<double(const std::vector<double>&)>& phi,
                                const std::vector<double>& x0,
                                const MinimizerOptions& opts = {});

}
