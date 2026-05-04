// Simple gradient-based minimizer utilities.
#ifndef PP_OPTIMIZER_H
#define PP_OPTIMIZER_H

#include <functional>
#include <vector>

#include "vector.h"

namespace pp {

// Controls the minimizer behavior.
struct MinimizeOptions {
  int max_iters = 2000;
  int max_line_search = 20;
  double grad_tol = 1e-6;
  double step_init = 0.05;
  double step_min = 1e-8;
  double step_shrink = 0.5;
  double step_grow = 1.1;
};

// Reports the outcome of a minimization run.
struct MinimizeReport {
  int iters = 0;
  double final_cost = 0.0;
  double final_grad_norm = 0.0;
  std::vector<double> cost_history;
};

// Cost function signature; pass nullptr to skip gradient evaluation.
using CostFunction = std::function<double(const pp::vector&, pp::vector*)>;

// Euclidean norm helper for gradients.
double l2_norm(const pp::vector& v);

// Minimizes the provided cost function in-place.
MinimizeReport minimize(pp::vector& params, const CostFunction& cost_fn,
                        const MinimizeOptions& opt,
                        std::vector<double>* cost_history = nullptr);

}  // namespace pp

#endif  // PP_OPTIMIZER_H
