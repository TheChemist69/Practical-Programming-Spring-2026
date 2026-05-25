// "ann.h" header file.
// One-input, one-output artificial neural network with analytic derivatives.
//
// Network form: F(x) = sum_{i=0}^{n-1} w_i * sigma( a_i * (x - b_i) )
// where sigma is the activation function, a_i are width parameters,
// b_i are center parameters, and w_i are output weights.
// For ODE training the second derivative F''(x) is also available analytically.
#ifndef PP_ANN_H
#define PP_ANN_H

#include "activation.h"
#include "optimizer.h"
#include "vector.h"

#include <functional>

namespace pp {

// Training options for tabulated function interpolation.
struct TrainOptions {
  double x_min = -1.0;            // left boundary of the training interval
  double x_max = 1.0;             // right boundary of the training interval
  double init_log_b = -0.7;       // log10 of the initial neuron width (center spread scale)
  unsigned seed = 1;              // random seed for parameter initialization
  bool use_coordinate_sweep = false; // if true, optimize one parameter at a time
  int max_sweeps = 50;            // maximum coordinate-sweep passes
  double sweep_tol = 1e-8;        // stop coordinate sweeps when improvement < sweep_tol
  MinimizeOptions minimize;       // options forwarded to the gradient minimizer
};

// Training options for differential-equation problems.
// The ODE cost is: integral of phi(F, F', F'')^2 dx + alpha*(F(c)-y_c)^2 + beta*(F'(c)-y1_c)^2
struct OdeOptions {
  double a = -1.0;               // left boundary of the ODE integration domain
  double b = 1.0;                // right boundary of the ODE integration domain
  double c = 0.0;                // point where boundary conditions are imposed
  double y_c = 0.0;              // required value F(c) = y_c
  double y1_c = 0.0;             // required value F'(c) = y1_c
  double alpha = 10.0;           // penalty weight for the F(c) = y_c condition
  double beta = 10.0;            // penalty weight for the F'(c) = y1_c condition
  int samples = 64;              // number of quadrature points for the ODE cost integral
  double init_log_b = -0.7;      // log10 of the initial neuron width
  unsigned seed = 1;             // random seed for parameter initialization
  bool use_coordinate_sweep = false; // if true, optimize one parameter at a time
  int max_sweeps = 50;           // maximum coordinate-sweep passes
  double sweep_tol = 1e-8;       // stop coordinate sweeps when improvement < sweep_tol
  MinimizeOptions minimize;      // options forwarded to the gradient minimizer
};

// Residual and partial derivatives of an ODE at a single evaluation point.
// phi(y, y', y'', x) = 0 describes the ODE; the partials are needed for the gradient.
struct OdeResidual {
  double phi = 0.0;      // ODE residual value at this point
  double dphi_dy = 0.0;  // partial d(phi)/d(y)
  double dphi_dy1 = 0.0; // partial d(phi)/d(y')
  double dphi_dy2 = 0.0; // partial d(phi)/d(y'')
};

// A function returning the ODE residual and its partial derivatives at (y, y', y'', x).
using OdeFunction = std::function<OdeResidual(double y, double y1, double y2, double x)>;

// Summary of a training run returned by Ann::train / Ann::train_ode.
struct TrainReport {
  MinimizeReport minimize; // iteration count, final cost, gradient norm, cost history
};

// One-input, one-output ANN: F(x) = sum_i w_i * sigma(a_i * (x - b_i))
class Ann {
 public:
  // Creates a network with n hidden neurons and the supplied activation function.
  Ann(int n, const Activation& activation);

  // Returns F(x).
  double response(double x) const;

  // Returns F'(x).
  double response_deriv(double x) const;

  // Returns F''(x).
  double response_second(double x) const;

  // Returns the definite integral from x0 to x of F.
  double response_antideriv(double x, double x0) const;

  // Initializes (a_i, b_i, w_i) with centers spread uniformly over [x_min, x_max].
  void initialize_uniform(double x_min, double x_max, double init_log_b, unsigned seed);

  // Trains the network to fit tabulated data (x, y) by minimizing sum of squared residuals.
  TrainReport train(const pp::vector& x, const pp::vector& y,
                    const TrainOptions& options);

  // Trains the network to satisfy an ODE on [a, b] subject to boundary conditions at c.
  TrainReport train_ode(const OdeFunction& ode, const OdeOptions& options);

  // Returns the raw parameter vector [a_0, log_b_0, w_0, a_1, log_b_1, w_1, ...].
  const pp::vector& parameters() const { return params_; }

 private:
  double compute_interpolation_cost(const pp::vector& params,
                                    const pp::vector& x, const pp::vector& y,
                                    pp::vector* grad) const;

  double compute_ode_cost(const pp::vector& params, const OdeFunction& ode,
                          const OdeOptions& options, pp::vector* grad) const;

  int n_ = 0;
  Activation activation_;
  pp::vector params_;  // layout: [a_0, log_b_0, w_0, a_1, log_b_1, w_1, ...]

  double a_at(int i) const;
  double log_b_at(int i) const;
  double b_at(int i) const;
  double w_at(int i) const;
};

}  // namespace pp

#endif  // PP_ANN_H
