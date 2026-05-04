// Artificial neural network interface and training helpers.
#ifndef PP_ANN_H
#define PP_ANN_H

#include "activation.h"
#include "optimizer.h"
#include "vector.h"

#include <functional>

namespace pp {

// Training options for tabulated interpolation.
struct TrainOptions {
  double x_min = -1.0;
  double x_max = 1.0;
  double init_log_b = -0.7;
  unsigned seed = 1;
  bool use_coordinate_sweep = false;
  int max_sweeps = 50;
  double sweep_tol = 1e-8;
  MinimizeOptions minimize;
};

// Training options for differential-equation problems.
struct OdeOptions {
  double a = -1.0;
  double b = 1.0;
  double c = 0.0;
  double y_c = 0.0;
  double y1_c = 0.0;
  double alpha = 10.0;
  double beta = 10.0;
  int samples = 64;
  double init_log_b = -0.7;
  unsigned seed = 1;
  bool use_coordinate_sweep = false;
  int max_sweeps = 50;
  double sweep_tol = 1e-8;
  MinimizeOptions minimize;
};

// Residual and partial derivatives for a differential equation.
struct OdeResidual {
  double phi = 0.0;
  double dphi_dy = 0.0;
  double dphi_dy1 = 0.0;
  double dphi_dy2 = 0.0;
};

using OdeFunction = std::function<OdeResidual(double y, double y1, double y2, double x)>;

// Wraps minimizer results for reporting.
struct TrainReport {
  MinimizeReport minimize;
};

// One-input, one-output ANN with shared activation and analytic derivatives.
class Ann {
 public:
  // Creates a network with n hidden neurons and the supplied activation.
  Ann(int n, const Activation& activation);

  // Returns the network response at x.
  double response(double x) const;

  // Returns the first derivative of the response at x.
  double response_deriv(double x) const;

  // Returns the second derivative of the response at x.
  double response_second(double x) const;

  // Returns the definite integral from x0 to x of the response.
  double response_antideriv(double x, double x0) const;

  // Initializes parameters using a uniform spread of centers.
  void initialize_uniform(double x_min, double x_max, double init_log_b, unsigned seed);

  // Trains the network on a tabulated function.
  TrainReport train(const pp::vector& x, const pp::vector& y,
                    const TrainOptions& options);

  // Trains the network to satisfy a differential equation on [a, b].
  TrainReport train_ode(const OdeFunction& ode, const OdeOptions& options);

  // Exposes the parameter vector for diagnostics.
  const pp::vector& parameters() const { return params_; }

 private:
  double compute_interpolation_cost(const pp::vector& params,
                                    const pp::vector& x, const pp::vector& y,
                                    pp::vector* grad) const;

  double compute_ode_cost(const pp::vector& params, const OdeFunction& ode,
                          const OdeOptions& options, pp::vector* grad) const;

  int n_ = 0;
  Activation activation_;
  pp::vector params_;

  double a_at(int i) const;
  double log_b_at(int i) const;
  double b_at(int i) const;
  double w_at(int i) const;
};

}  // namespace pp

#endif  // PP_ANN_H
