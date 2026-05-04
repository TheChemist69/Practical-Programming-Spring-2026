// Gradient-based minimizer implementation.
#include "optimizer.h"

#include <algorithm>
#include <cmath>

namespace pp {

// Returns the Euclidean norm of a vector.
double l2_norm(const pp::vector& v) {
  double sum = 0.0;
  for (int i = 0; i < v.size(); ++i) {
    const double value = v[i];
    sum += value * value;
  }
  return std::sqrt(sum);
}

// Runs a simple backtracking line-search gradient descent.
MinimizeReport minimize(pp::vector& params, const CostFunction& cost_fn,
                        const MinimizeOptions& opt,
                        std::vector<double>* cost_history) {
  pp::vector grad(params.size());
  pp::vector trial(params.size());

  double cost = cost_fn(params, &grad);
  if (cost_history) {
    cost_history->clear();
    cost_history->push_back(cost);
  }
  double grad_norm = l2_norm(grad);
  double step = opt.step_init;

  int iter = 0;
  for (; iter < opt.max_iters; ++iter) {
    if (grad_norm <= opt.grad_tol) {
      break;
    }

    bool accepted = false;
    double local_step = step;
    double best_cost = cost;

    for (int ls = 0; ls < opt.max_line_search; ++ls) {
      for (int i = 0; i < params.size(); ++i) {
        trial[i] = params[i] - local_step * grad[i];
      }

      const double trial_cost = cost_fn(trial, nullptr);
      if (trial_cost < best_cost) {
        params = trial;
        best_cost = trial_cost;
        accepted = true;
        step = std::max(opt.step_min, local_step * opt.step_grow);
        break;
      }

      local_step *= opt.step_shrink;
      if (local_step < opt.step_min) {
        break;
      }
    }

    if (!accepted) {
      step = std::max(opt.step_min, step * opt.step_shrink);
      if (step <= opt.step_min) {
        break;
      }
    }

    cost = cost_fn(params, &grad);
    if (cost_history) {
      cost_history->push_back(cost);
    }
    grad_norm = l2_norm(grad);
  }

  MinimizeReport report;
  if (cost_history) {
    report.cost_history = *cost_history;
  }
  report.iters = iter;
  report.final_cost = cost;
  report.final_grad_norm = grad_norm;
  return report;
}

}  // namespace pp
