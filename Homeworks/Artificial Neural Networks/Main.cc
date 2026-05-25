// Main driver for the ANN homework.
#include "activation.h"
#include "ann.h"
#include "output_writer.h"
#include "vector.h"

#include <cstddef>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace {

// Returns the target function g(x) = cos(5x - 1) * exp(-x^2).
double target_function(double x) {
  return std::cos(5.0 * x - 1.0) * std::exp(-x * x);
}

// Returns the first derivative of g(x).
double target_function_d1(double x) {
  const double phase = 5.0 * x - 1.0;
  const double e = std::exp(-x * x);
  const double s = std::sin(phase);
  const double c = std::cos(phase);
  return e * (-5.0 * s - 2.0 * x * c);
}

// Returns the second derivative of g(x).
double target_function_d2(double x) {
  const double phase = 5.0 * x - 1.0;
  const double e = std::exp(-x * x);
  const double s = std::sin(phase);
  const double c = std::cos(phase);
  return e * (20.0 * x * s + (4.0 * x * x - 27.0) * c);
}

// Creates an evenly spaced grid of n points from a to b.
pp::vector linspace(double a, double b, int n) {
  pp::vector xs;
  if (n <= 0) {
    return xs;
  }
  xs.resize(n);
  if (n == 1) {
    xs[0] = a;
    return xs;
  }

  const double step = (b - a) / static_cast<double>(n - 1);
  for (int i = 0; i < n; ++i) {
    xs[i] = a + step * static_cast<double>(i);
  }
  return xs;
}

// Computes cumulative trapezoid integral from x[0] to each x[i].
pp::vector cumulative_trapezoid(const pp::vector& x, const pp::vector& y) {
  const int n = x.size();
  pp::vector result(n);
  if (n < 2) {
    return result;
  }

  for (int i = 1; i < n; ++i) {
    const double dx = x[i] - x[i - 1];
    result[i] = result[i - 1] + 0.5 * dx * (y[i] + y[i - 1]);
  }
  return result;
}

// Computes the RMS error between two vectors.
double rms_error(const pp::vector& a, const pp::vector& b) {
  if (a.size() != b.size() || a.size() == 0) {
    return 0.0;
  }
  double sum = 0.0;
  for (int i = 0; i < a.size(); ++i) {
    const double diff = a[i] - b[i];
    sum += diff * diff;
  }
  return std::sqrt(sum / static_cast<double>(a.size()));
}

}  // namespace

// Runs the three exercises and writes Out.txtx.
int main() {
  const double x_min = -1.0;
  const double x_max = 1.0;
  const int neurons = 8;
  const int samples = 40;

  std::ofstream out("Out.txtx");
  if (!out) {
    std::cerr << "Failed to open Out.txtx for writing.\n";
    return 1;
  }

  std::cout << "Running ANN homework. Writing Out.txtx in this folder." << std::endl;

  out << std::scientific << std::setprecision(8);
  const int col_w = 16;

  report::write_title(out, "Practical Programming - Homework: Artificial Neural Networks");
  out << "# Network: Fp(x) = sum_i w_i f((x-a_i)/b_i)\n";
  out << "# Cost: C(p) = sum_k (Fp(x_k) - y_k)^2\n";
  out << "# Sections: (1) interpolation, (2) derivatives/antiderivative, (3) ODE training\n";
  report::write_rule(out);
  out << "\n";

  // Exercise 1: interpolation of the tabulated function.
  const pp::Activation act = pp::gaussian_wavelet_activation();
  pp::Ann net(neurons, act);

  const pp::vector xs = linspace(x_min, x_max, samples);
  pp::vector ys(xs.size());
  for (int i = 0; i < xs.size(); ++i) {
    ys[i] = target_function(xs[i]);
  }

  pp::TrainOptions train_opts;
  train_opts.x_min = x_min;
  train_opts.x_max = x_max;
  train_opts.init_log_b = std::log((x_max - x_min) / static_cast<double>(neurons));
  train_opts.seed = 2;
  train_opts.use_coordinate_sweep = true;
  train_opts.max_sweeps = 40;
  train_opts.sweep_tol = 1e-8;
  train_opts.minimize.max_iters = 4000;
  train_opts.minimize.grad_tol = 1e-6;
  train_opts.minimize.step_init = 0.05;
  train_opts.minimize.step_min = 1e-8;

  const pp::TrainReport report = net.train(xs, ys, train_opts);

  pp::vector y_fit(xs.size());
  for (int i = 0; i < xs.size(); ++i) {
    y_fit[i] = net.response(xs[i]);
  }

  const double rms_fit = rms_error(y_fit, ys);

  std::cout << "Exercise 1 done. RMS error = " << rms_fit << std::endl;

  std::cout << "Writing ann_train_history.data" << std::endl;
  {
    std::ofstream history_data("ann_train_history.data");
    if (!history_data) {
      std::cerr << "Failed to open ann_train_history.data for writing.\n";
      return 1;
    }
    history_data << std::scientific << std::setprecision(8);
    history_data << "# iter cost\n";
    for (std::size_t i = 0; i < report.minimize.cost_history.size(); ++i) {
      history_data << i << " " << report.minimize.cost_history[i] << "\n";
    }
  }

  std::cout << "Writing ann_fit.data" << std::endl;
  {
    std::ofstream fit_data("ann_fit.data");
    if (!fit_data) {
      std::cerr << "Failed to open ann_fit.data for writing.\n";
      return 1;
    }

    fit_data << std::scientific << std::setprecision(8);
    fit_data << "# x target ann\n";
    for (int i = 0; i < xs.size(); ++i) {
      fit_data << xs[i] << " " << ys[i] << " " << y_fit[i] << "\n";
    }
  }

  std::cout << "Writing ann_error.data" << std::endl;
  {
    std::ofstream error_data("ann_error.data");
    if (!error_data) {
      std::cerr << "Failed to open ann_error.data for writing.\n";
      return 1;
    }

    error_data << std::scientific << std::setprecision(8);
    error_data << "# x target ann error abs_error\n";
    for (int i = 0; i < xs.size(); ++i) {
      const double error = y_fit[i] - ys[i];
      error_data << xs[i] << " " << ys[i] << " " << y_fit[i] << " " << error
                 << " " << std::abs(error) << "\n";
    }
  }

  report::write_title(out, "Exercise 1 - Interpolation of g(x)");
  report::write_kv(out, "Target", "g(x)=cos(5x-1)exp(-x^2)");
  report::write_kv(out, "Activation", act.name.c_str());
  report::write_kv(out, "Hidden neurons", neurons);
  report::write_interval(out, "Sample range", x_min, x_max);
  report::write_kv(out, "Samples", samples);
  report::write_kv(out, "Init log(b)", train_opts.init_log_b);
  report::write_kv(out, "Max iters", train_opts.minimize.max_iters);
  report::write_kv(out, "Grad tol", train_opts.minimize.grad_tol);
  report::write_kv(out, "Step init", train_opts.minimize.step_init);
  report::write_kv(out, "Step min", train_opts.minimize.step_min);
  report::write_kv(out, "Iterations", report.minimize.iters);
  report::write_kv(out, "Final cost", report.minimize.final_cost);
  report::write_kv(out, "Final grad norm", report.minimize.final_grad_norm);
  report::write_kv(out, "RMS error", rms_fit);
  out << "#\n";
  out << "# " << std::left << std::setw(col_w) << "x" << std::setw(col_w) << "g(x)"
      << std::setw(col_w) << "Fp(x)" << "\n";
  out << std::right;
  for (int i = 0; i < xs.size(); ++i) {
    out << std::setw(col_w) << xs[i] << std::setw(col_w) << ys[i]
        << std::setw(col_w) << y_fit[i] << "\n";
  }
  out << "\n";

  // Exercise 2: derivatives and antiderivative of the approximant.
  pp::vector g_d1(xs.size());
  pp::vector g_d2(xs.size());
  pp::vector ann_d1(xs.size());
  pp::vector ann_d2(xs.size());
  pp::vector ann_int(xs.size());

  for (int i = 0; i < xs.size(); ++i) {
    const double x = xs[i];
    g_d1[i] = target_function_d1(x);
    g_d2[i] = target_function_d2(x);
    ann_d1[i] = net.response_deriv(x);
    ann_d2[i] = net.response_second(x);
    ann_int[i] = net.response_antideriv(x, xs[0]);
  }

  const pp::vector g_int = cumulative_trapezoid(xs, ys);
  const double rms_d1 = rms_error(ann_d1, g_d1);
  const double rms_d2 = rms_error(ann_d2, g_d2);

  std::cout << "Exercise 2 done. RMS g' = " << rms_d1
            << ", RMS g'' = " << rms_d2 << std::endl;

  std::cout << "Writing ann_derivatives.data" << std::endl;
  {
    std::ofstream deriv_data("ann_derivatives.data");
    if (!deriv_data) {
      std::cerr << "Failed to open ann_derivatives.data for writing.\n";
      return 1;
    }

    deriv_data << std::scientific << std::setprecision(8);
    deriv_data << "# x g1 ann1 g2 ann2 gint annint\n";
    for (int i = 0; i < xs.size(); ++i) {
      deriv_data << xs[i] << " " << g_d1[i] << " " << ann_d1[i] << " "
                 << g_d2[i] << " " << ann_d2[i] << " " << g_int[i] << " "
                 << ann_int[i] << "\n";
    }
  }

  report::write_title(out, "Exercise 2 - Derivatives and Antiderivative");
  report::write_kv(out, "Antiderivative", "Integral from x0 to x, x0 = first sample");
  report::write_kv(out, "RMS g'(x)", rms_d1);
  report::write_kv(out, "RMS g''(x)", rms_d2);
  out << "#\n";
  out << "# " << std::left << std::setw(col_w) << "x" << std::setw(col_w)
      << "g'(x)" << std::setw(col_w) << "Fp'(x)" << std::setw(col_w)
      << "g''(x)" << std::setw(col_w) << "Fp''(x)" << std::setw(col_w)
      << "Int g" << std::setw(col_w) << "Int Fp" << "\n";
  out << std::right;
  for (int i = 0; i < xs.size(); ++i) {
    out << std::setw(col_w) << xs[i] << std::setw(col_w) << g_d1[i]
        << std::setw(col_w) << ann_d1[i] << std::setw(col_w) << g_d2[i]
        << std::setw(col_w) << ann_d2[i] << std::setw(col_w) << g_int[i]
        << std::setw(col_w) << ann_int[i] << "\n";
  }
  out << "\n";

  // Exercise 3: ODE y'' + y = 0 with y(0)=0 and y'(0)=1.
  pp::Ann ode_net(neurons, act);
  pp::OdeOptions ode_opts;
  ode_opts.a = -1.0;
  ode_opts.b = 1.0;
  ode_opts.c = 0.0;
  ode_opts.y_c = 0.0;
  ode_opts.y1_c = 1.0;
  ode_opts.alpha = 10.0;
  ode_opts.beta = 10.0;
  ode_opts.samples = 80;
  ode_opts.use_coordinate_sweep = true;
  ode_opts.max_sweeps = 40;
  ode_opts.sweep_tol = 1e-8;
  ode_opts.init_log_b = std::log((ode_opts.b - ode_opts.a) / static_cast<double>(neurons));
  ode_opts.seed = 3;
  ode_opts.minimize.max_iters = 5000;
  ode_opts.minimize.grad_tol = 1e-6;
  ode_opts.minimize.step_init = 0.03;
  ode_opts.minimize.step_min = 1e-8;

  const pp::OdeFunction ode = [](double y, double /*y1*/, double y2, double /*x*/) {
    pp::OdeResidual res;
    res.phi = y2 + y;
    res.dphi_dy = 1.0;
    res.dphi_dy1 = 0.0;
    res.dphi_dy2 = 1.0;
    return res;
  };

  const pp::TrainReport ode_report = ode_net.train_ode(ode, ode_opts);

  std::cout << "Writing ann_ode_history.data" << std::endl;
  {
    std::ofstream ode_history("ann_ode_history.data");
    if (!ode_history) {
      std::cerr << "Failed to open ann_ode_history.data for writing.\n";
      return 1;
    }
    ode_history << std::scientific << std::setprecision(8);
    ode_history << "# iter cost\n";
    for (std::size_t i = 0; i < ode_report.minimize.cost_history.size(); ++i) {
      ode_history << i << " " << ode_report.minimize.cost_history[i] << "\n";
    }
  }

  const pp::vector xs_ode = linspace(ode_opts.a, ode_opts.b, 60);
  pp::vector y_true(xs_ode.size());
  pp::vector y_pred(xs_ode.size());
  pp::vector y_pred_d1(xs_ode.size());

  for (int i = 0; i < xs_ode.size(); ++i) {
    const double x = xs_ode[i];
    y_true[i] = std::sin(x);
    y_pred[i] = ode_net.response(x);
    y_pred_d1[i] = ode_net.response_deriv(x);
  }

  const double rms_ode = rms_error(y_pred, y_true);

  std::cout << "Exercise 3 done. RMS vs sin(x) = " << rms_ode << std::endl;

  std::cout << "Writing ann_ode.data" << std::endl;
  {
    std::ofstream ode_data("ann_ode.data");
    if (!ode_data) {
      std::cerr << "Failed to open ann_ode.data for writing.\n";
      return 1;
    }

    ode_data << std::scientific << std::setprecision(8);
    ode_data << "# x true ann ann_d1 true_d1\n";
    for (int i = 0; i < xs_ode.size(); ++i) {
      const double true_d1 = std::cos(xs_ode[i]);
      ode_data << xs_ode[i] << " " << y_true[i] << " " << y_pred[i] << " "
               << y_pred_d1[i] << " " << true_d1 << "\n";
    }
  }

  std::cout << "Writing ann_ode_residual.data" << std::endl;
  {
    std::ofstream ode_residual("ann_ode_residual.data");
    if (!ode_residual) {
      std::cerr << "Failed to open ann_ode_residual.data for writing.\n";
      return 1;
    }

    ode_residual << std::scientific << std::setprecision(8);
    ode_residual << "# x phi phi2 abs_error\n";
    for (int i = 0; i < xs_ode.size(); ++i) {
      const double x = xs_ode[i];
      const double y2 = ode_net.response_second(x);
      const pp::OdeResidual res = ode(y_pred[i], y_pred_d1[i], y2, x);
      const double phi = res.phi;
      const double err = y_pred[i] - y_true[i];
      ode_residual << x << " " << phi << " " << phi * phi << " " << std::abs(err)
                   << "\n";
    }
  }

  report::write_title(out, "Exercise 3 - ODE Training");
  report::write_kv(out, "Equation", "y'' + y = 0");
  report::write_kv(out, "Boundary", "y(0)=0, y'(0)=1");
  report::write_interval(out, "Interval", ode_opts.a, ode_opts.b);
  report::write_kv(out, "c", ode_opts.c);
  report::write_kv(out, "alpha", ode_opts.alpha);
  report::write_kv(out, "beta", ode_opts.beta);
  report::write_kv(out, "Samples", ode_opts.samples);
  report::write_kv(out, "Activation", act.name.c_str());
  report::write_kv(out, "Hidden neurons", neurons);
  report::write_kv(out, "Iterations", ode_report.minimize.iters);
  report::write_kv(out, "Final cost", ode_report.minimize.final_cost);
  report::write_kv(out, "Final grad norm", ode_report.minimize.final_grad_norm);
  report::write_kv(out, "RMS error vs sin(x)", rms_ode);
  out << "#\n";
  out << "# Cost: integral_a^b Phi^2 dx + alpha(Fp(c)-yc)^2 + beta(Fp'(c)-y'c)^2\n";
  out << "# Phi(y'',y',y,x) = y'' + y\n";
  out << "#\n";
  out << "# " << std::left << std::setw(col_w) << "x" << std::setw(col_w)
      << "y_true" << std::setw(col_w) << "Fp(x)" << std::setw(col_w)
      << "Fp'(x)" << "\n";
  out << std::right;
  for (int i = 0; i < xs_ode.size(); ++i) {
    out << std::setw(col_w) << xs_ode[i] << std::setw(col_w) << y_true[i]
        << std::setw(col_w) << y_pred[i] << std::setw(col_w) << y_pred_d1[i]
        << "\n";
  }

  out.flush();
  std::cout << "Finished. Out.txtx is ready." << std::endl;

  return 0;
}
