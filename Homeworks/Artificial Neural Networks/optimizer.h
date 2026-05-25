// "optimizer.h" header file.
// Simple gradient-descent minimizer with adaptive step-size (Barzilai-Borwein style).
#ifndef PP_OPTIMIZER_H
#define PP_OPTIMIZER_H

#include <functional>
#include <vector>

#include "vector.h"

namespace pp {

// Controls the gradient-descent minimizer behavior.
struct MinimizeOptions {
  int max_iters = 2000;       // maximum number of gradient steps
  int max_line_search = 20;   // maximum number of backtracking halvings per step
  double grad_tol = 1e-6;     // convergence criterion: stop when ||grad|| < grad_tol
  double step_init = 0.05;    // initial step size
  double step_min = 1e-8;     // minimum allowed step size before giving up
  double step_shrink = 0.5;   // multiplicative factor for step reduction during backtracking
  double step_grow = 1.1;     // multiplicative factor for step growth after a successful step
};

// Reports the outcome of a minimization run.
struct MinimizeReport {
  int iters = 0;                        // number of gradient steps taken
  double final_cost = 0.0;             // cost at the returned parameters
  double final_grad_norm = 0.0;        // ||gradient|| at the returned parameters
  std::vector<double> cost_history;    // cost value after each iteration (for plotting convergence)
};

// Cost function signature: evaluates cost and optionally writes the gradient.
// If grad != nullptr, the function must fill *grad; otherwise it may skip it.
using CostFunction = std::function<double(const pp::vector&, pp::vector*)>;

// Euclidean norm of a vector (used internally for gradient convergence checks).
double l2_norm(const pp::vector& v);

// Runs gradient descent on params in-place, using cost_fn for gradient information.
// Optionally appends cost values to cost_history if non-null.
MinimizeReport minimize(pp::vector& params, const CostFunction& cost_fn,
                        const MinimizeOptions& opt,
                        std::vector<double>* cost_history = nullptr);

}  // namespace pp

#endif  // PP_OPTIMIZER_H
