// Artificial neural network implementation.
#include "ann.h"

#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace pp {
namespace {

// Computes a trapezoid integration weight for index i in [0, n-1].
double trapezoid_weight(int i, int n, double dx) {
  if (n < 2) {
    return 0.0;
  }
  if (i == 0 || i == n - 1) {
    return 0.5 * dx;
  }
  return dx;
}

}  // namespace

// Constructs the network with n hidden neurons.
Ann::Ann(int n, const Activation& activation)
  : n_(n), activation_(activation), params_(3 * n) {
  if (n_ <= 0) {
    throw std::invalid_argument("Ann: n must be positive");
  }
}

// Returns the response at x.
double Ann::response(double x) const {
  double sum = 0.0;
  for (int i = 0; i < n_; ++i) {
    const double a = a_at(i);
    const double b = b_at(i);
    const double w = w_at(i);
    const double z = (x - a) / b;
    sum += w * activation_.f(z);
  }
  return sum;
}

// Returns the first derivative at x.
double Ann::response_deriv(double x) const {
  double sum = 0.0;
  for (int i = 0; i < n_; ++i) {
    const double a = a_at(i);
    const double b = b_at(i);
    const double w = w_at(i);
    const double z = (x - a) / b;
    sum += w * activation_.df(z) / b;
  }
  return sum;
}

// Returns the second derivative at x.
double Ann::response_second(double x) const {
  double sum = 0.0;
  for (int i = 0; i < n_; ++i) {
    const double a = a_at(i);
    const double b = b_at(i);
    const double w = w_at(i);
    const double z = (x - a) / b;
    sum += w * activation_.ddf(z) / (b * b);
  }
  return sum;
}

// Returns the definite integral from x0 to x.
double Ann::response_antideriv(double x, double x0) const {
  double sum = 0.0;
  for (int i = 0; i < n_; ++i) {
    const double a = a_at(i);
    const double b = b_at(i);
    const double w = w_at(i);
    const double z = (x - a) / b;
    const double z0 = (x0 - a) / b;
    sum += w * b * (activation_.antideriv(z) - activation_.antideriv(z0));
  }
  return sum;
}

// Initializes parameters with evenly spaced centers.
void Ann::initialize_uniform(double x_min, double x_max, double init_log_b,
                             unsigned seed) {
  const double span = x_max - x_min;
  const double step = span / static_cast<double>(n_);

  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  for (int i = 0; i < n_; ++i) {
    params_[3 * i] = x_min + (i + 0.5) * step;
    params_[3 * i + 1] = init_log_b;
    params_[3 * i + 2] = dist(rng);
  }
}

// Computes the interpolation cost and optional gradient for the given params.
double Ann::compute_interpolation_cost(const pp::vector& params,
                                       const pp::vector& x, const pp::vector& y,
                                       pp::vector* grad) const {
  double cost = 0.0;
  if (grad) {
    grad->resize(params.size());
    for (int i = 0; i < grad->size(); ++i) {
      (*grad)[i] = 0.0;
    }
  }

  for (int k = 0; k < x.size(); ++k) {
    double sum = 0.0;
    for (int i = 0; i < n_; ++i) {
      const double a = params[3 * i];
      const double log_b = params[3 * i + 1];
      const double b = std::exp(log_b);
      const double w = params[3 * i + 2];
      const double z = (x[k] - a) / b;
      sum += w * activation_.f(z);
    }

    const double residual = sum - y[k];
    cost += residual * residual;

    if (!grad) {
      continue;
    }

    for (int i = 0; i < n_; ++i) {
      const double a = params[3 * i];
      const double log_b = params[3 * i + 1];
      const double b = std::exp(log_b);
      const double w = params[3 * i + 2];
      const double z = (x[k] - a) / b;
      const double f0 = activation_.f(z);
      const double f1 = activation_.df(z);

      const double d_f_da = -w * f1 / b;
      const double d_f_dlogb = -w * f1 * z;
      const double d_f_dw = f0;

      (*grad)[3 * i] += 2.0 * residual * d_f_da;
      (*grad)[3 * i + 1] += 2.0 * residual * d_f_dlogb;
      (*grad)[3 * i + 2] += 2.0 * residual * d_f_dw;
    }
  }

  return cost;
}

// Computes the ODE cost and optional gradient for the given params.
double Ann::compute_ode_cost(const pp::vector& params, const OdeFunction& ode,
                             const OdeOptions& options, pp::vector* grad) const {
  const double a = options.a;
  const double b = options.b;
  const double dx = (b - a) / static_cast<double>(options.samples - 1);

  double cost = 0.0;
  if (grad) {
    grad->resize(params.size());
    for (int i = 0; i < grad->size(); ++i) {
      (*grad)[i] = 0.0;
    }
  }

  for (int k = 0; k < options.samples; ++k) {
    const double x = a + dx * static_cast<double>(k);
    double y = 0.0;
    double y1 = 0.0;
    double y2 = 0.0;

    for (int i = 0; i < n_; ++i) {
      const double a_i = params[3 * i];
      const double log_b = params[3 * i + 1];
      const double b_i = std::exp(log_b);
      const double w_i = params[3 * i + 2];
      const double z = (x - a_i) / b_i;

      y += w_i * activation_.f(z);
      y1 += w_i * activation_.df(z) / b_i;
      y2 += w_i * activation_.ddf(z) / (b_i * b_i);
    }

    const OdeResidual res = ode(y, y1, y2, x);
    const double weight = trapezoid_weight(k, options.samples, dx);
    cost += weight * res.phi * res.phi;

    if (!grad) {
      continue;
    }

    for (int i = 0; i < n_; ++i) {
      const double a_i = params[3 * i];
      const double log_b = params[3 * i + 1];
      const double b_i = std::exp(log_b);
      const double w_i = params[3 * i + 2];
      const double z = (x - a_i) / b_i;

      const double f0 = activation_.f(z);
      const double f1 = activation_.df(z);
      const double f2 = activation_.ddf(z);
      const double f3 = activation_.dddf(z);

      const double d_f_da = -w_i * f1 / b_i;
      const double d_f_dlogb = -w_i * f1 * z;
      const double d_f_dw = f0;

      const double d_f1_da = -w_i * f2 / (b_i * b_i);
      const double d_f1_dlogb = -w_i * (f1 + z * f2) / b_i;
      const double d_f1_dw = f1 / b_i;

      const double d_f2_da = -w_i * f3 / (b_i * b_i * b_i);
      const double d_f2_dlogb = -w_i * (2.0 * f2 + z * f3) / (b_i * b_i);
      const double d_f2_dw = f2 / (b_i * b_i);

      const double d_phi = res.dphi_dy * d_f_da + res.dphi_dy1 * d_f1_da +
                           res.dphi_dy2 * d_f2_da;
      const double d_phi_logb =
          res.dphi_dy * d_f_dlogb + res.dphi_dy1 * d_f1_dlogb +
          res.dphi_dy2 * d_f2_dlogb;
      const double d_phi_w =
          res.dphi_dy * d_f_dw + res.dphi_dy1 * d_f1_dw + res.dphi_dy2 * d_f2_dw;

      (*grad)[3 * i] += 2.0 * weight * res.phi * d_phi;
      (*grad)[3 * i + 1] += 2.0 * weight * res.phi * d_phi_logb;
      (*grad)[3 * i + 2] += 2.0 * weight * res.phi * d_phi_w;
    }
  }

  double yc = 0.0;
  double y1c = 0.0;
  for (int i = 0; i < n_; ++i) {
    const double a_i = params[3 * i];
    const double log_b = params[3 * i + 1];
    const double b_i = std::exp(log_b);
    const double w_i = params[3 * i + 2];
    const double z = (options.c - a_i) / b_i;

    yc += w_i * activation_.f(z);
    y1c += w_i * activation_.df(z) / b_i;
  }

  const double bc_res = yc - options.y_c;
  const double bc_res1 = y1c - options.y1_c;
  cost += options.alpha * bc_res * bc_res + options.beta * bc_res1 * bc_res1;

  if (grad) {
    for (int i = 0; i < n_; ++i) {
      const double a_i = params[3 * i];
      const double log_b = params[3 * i + 1];
      const double b_i = std::exp(log_b);
      const double w_i = params[3 * i + 2];
      const double z = (options.c - a_i) / b_i;

      const double f0 = activation_.f(z);
      const double f1 = activation_.df(z);
      const double f2 = activation_.ddf(z);

      const double d_f_da = -w_i * f1 / b_i;
      const double d_f_dlogb = -w_i * f1 * z;
      const double d_f_dw = f0;

      const double d_f1_da = -w_i * f2 / (b_i * b_i);
      const double d_f1_dlogb = -w_i * (f1 + z * f2) / b_i;
      const double d_f1_dw = f1 / b_i;

      (*grad)[3 * i] += 2.0 * options.alpha * bc_res * d_f_da +
                        2.0 * options.beta * bc_res1 * d_f1_da;
      (*grad)[3 * i + 1] += 2.0 * options.alpha * bc_res * d_f_dlogb +
                            2.0 * options.beta * bc_res1 * d_f1_dlogb;
      (*grad)[3 * i + 2] += 2.0 * options.alpha * bc_res * d_f_dw +
                            2.0 * options.beta * bc_res1 * d_f1_dw;
    }
  }

  return cost;
}

// Trains the network on tabulated data.
TrainReport Ann::train(const pp::vector& x, const pp::vector& y,
                       const TrainOptions& options) {
  if (x.size() != y.size()) {
    throw std::invalid_argument("Ann::train: x and y sizes must match");
  }

  initialize_uniform(options.x_min, options.x_max, options.init_log_b, options.seed);

  TrainReport result;
  if (!options.use_coordinate_sweep) {
    const CostFunction cost_fn = [&](const pp::vector& params, pp::vector* grad) {
      return compute_interpolation_cost(params, x, y, grad);
    };

    std::vector<double> cost_history;
    result.minimize = minimize(params_, cost_fn, options.minimize, &cost_history);
    return result;
  }

  std::vector<double> cost_history;
  double prev_cost = std::numeric_limits<double>::infinity();
  int sweep = 0;
  for (; sweep < options.max_sweeps; ++sweep) {
    for (int neuron = 0; neuron < n_; ++neuron) {
      pp::vector local(3);
      local[0] = params_[3 * neuron];
      local[1] = params_[3 * neuron + 1];
      local[2] = params_[3 * neuron + 2];

      const CostFunction local_cost = [&](const pp::vector& local_params,
                                          pp::vector* local_grad) {
        pp::vector trial = params_;
        trial[3 * neuron] = local_params[0];
        trial[3 * neuron + 1] = local_params[1];
        trial[3 * neuron + 2] = local_params[2];

        pp::vector full_grad;
        const double cost = compute_interpolation_cost(trial, x, y, &full_grad);
        if (local_grad) {
          local_grad->resize(3);
          (*local_grad)[0] = full_grad[3 * neuron];
          (*local_grad)[1] = full_grad[3 * neuron + 1];
          (*local_grad)[2] = full_grad[3 * neuron + 2];
        }
        return cost;
      };

      minimize(local, local_cost, options.minimize);
      params_[3 * neuron] = local[0];
      params_[3 * neuron + 1] = local[1];
      params_[3 * neuron + 2] = local[2];
    }

    const double cost = compute_interpolation_cost(params_, x, y, nullptr);
    cost_history.push_back(cost);
    if (std::abs(prev_cost - cost) <= options.sweep_tol) {
      break;
    }
    prev_cost = cost;
  }

  pp::vector final_grad;
  result.minimize.iters = sweep;
  result.minimize.cost_history = cost_history;
  result.minimize.final_cost = compute_interpolation_cost(params_, x, y, &final_grad);
  result.minimize.final_grad_norm = l2_norm(final_grad);
  return result;
}

// Trains the network to satisfy a differential equation.
TrainReport Ann::train_ode(const OdeFunction& ode, const OdeOptions& options) {
  if (options.samples < 2) {
    throw std::invalid_argument("Ann::train_ode: samples must be >= 2");
  }

  initialize_uniform(options.a, options.b, options.init_log_b, options.seed);

  TrainReport result;
  if (!options.use_coordinate_sweep) {
    const CostFunction cost_fn = [&](const pp::vector& params, pp::vector* grad) {
      return compute_ode_cost(params, ode, options, grad);
    };

    std::vector<double> cost_history;
    result.minimize = minimize(params_, cost_fn, options.minimize, &cost_history);
    return result;
  }

  std::vector<double> cost_history;
  double prev_cost = std::numeric_limits<double>::infinity();
  int sweep = 0;
  for (; sweep < options.max_sweeps; ++sweep) {
    for (int neuron = 0; neuron < n_; ++neuron) {
      pp::vector local(3);
      local[0] = params_[3 * neuron];
      local[1] = params_[3 * neuron + 1];
      local[2] = params_[3 * neuron + 2];

      const CostFunction local_cost = [&](const pp::vector& local_params,
                                          pp::vector* local_grad) {
        pp::vector trial = params_;
        trial[3 * neuron] = local_params[0];
        trial[3 * neuron + 1] = local_params[1];
        trial[3 * neuron + 2] = local_params[2];

        pp::vector full_grad;
        const double cost = compute_ode_cost(trial, ode, options, &full_grad);
        if (local_grad) {
          local_grad->resize(3);
          (*local_grad)[0] = full_grad[3 * neuron];
          (*local_grad)[1] = full_grad[3 * neuron + 1];
          (*local_grad)[2] = full_grad[3 * neuron + 2];
        }
        return cost;
      };

      minimize(local, local_cost, options.minimize);
      params_[3 * neuron] = local[0];
      params_[3 * neuron + 1] = local[1];
      params_[3 * neuron + 2] = local[2];
    }

    const double cost = compute_ode_cost(params_, ode, options, nullptr);
    cost_history.push_back(cost);
    if (std::abs(prev_cost - cost) <= options.sweep_tol) {
      break;
    }
    prev_cost = cost;
  }

  pp::vector final_grad;
  result.minimize.iters = sweep;
  result.minimize.cost_history = cost_history;
  result.minimize.final_cost = compute_ode_cost(params_, ode, options, &final_grad);
  result.minimize.final_grad_norm = l2_norm(final_grad);
  return result;
}
// Returns the i-th neuron center.
double Ann::a_at(int i) const {
  return params_[3 * i];
}

// Returns the i-th neuron log(b).
double Ann::log_b_at(int i) const {
  return params_[3 * i + 1];
}

// Returns the i-th neuron width.
double Ann::b_at(int i) const {
  return std::exp(log_b_at(i));
}

// Returns the i-th neuron weight.
double Ann::w_at(int i) const {
  return params_[3 * i + 2];
}

}  // namespace pp
