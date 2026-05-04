// Activation function implementations.
#include "activation.h"

#include <cmath>

namespace pp {
namespace {

// Gaussian wavelet and derivatives.
double gaussian_wavelet(double x) {
  const double x2 = x * x;
  return x * std::exp(-x2);
}

double gaussian_wavelet_d1(double x) {
  const double x2 = x * x;
  return std::exp(-x2) * (1.0 - 2.0 * x2);
}

double gaussian_wavelet_d2(double x) {
  const double x2 = x * x;
  const double x3 = x2 * x;
  return std::exp(-x2) * (4.0 * x3 - 6.0 * x);
}

double gaussian_wavelet_d3(double x) {
  const double x2 = x * x;
  const double x4 = x2 * x2;
  return std::exp(-x2) * (-8.0 * x4 + 24.0 * x2 - 6.0);
}

double gaussian_wavelet_antideriv(double x) {
  const double x2 = x * x;
  return -0.5 * std::exp(-x2);
}

// Gaussian and derivatives.
double gaussian(double x) {
  const double x2 = x * x;
  return std::exp(-x2);
}

double gaussian_d1(double x) {
  const double x2 = x * x;
  return -2.0 * x * std::exp(-x2);
}

double gaussian_d2(double x) {
  const double x2 = x * x;
  return (4.0 * x2 - 2.0) * std::exp(-x2);
}

double gaussian_d3(double x) {
  const double x2 = x * x;
  return (-8.0 * x2 * x + 12.0 * x) * std::exp(-x2);
}

double gaussian_antideriv(double x) {
  const double pi = std::acos(-1.0);
  return 0.5 * std::sqrt(pi) * std::erf(x);
}

}  // namespace

Activation gaussian_wavelet_activation() {
  Activation act;
  act.name = "gaussian_wavelet";
  act.f = gaussian_wavelet;
  act.df = gaussian_wavelet_d1;
  act.ddf = gaussian_wavelet_d2;
  act.dddf = gaussian_wavelet_d3;
  act.antideriv = gaussian_wavelet_antideriv;
  return act;
}

Activation gaussian_activation() {
  Activation act;
  act.name = "gaussian";
  act.f = gaussian;
  act.df = gaussian_d1;
  act.ddf = gaussian_d2;
  act.dddf = gaussian_d3;
  act.antideriv = gaussian_antideriv;
  return act;
}

}  // namespace pp
