// "activation.h" header file.
// Activation function bundles for the ANN.
//
// Each Activation carries the function value and its first three derivatives
// plus its antiderivative, all as std::function objects.  Having all four
// derivatives pre-computed avoids finite-difference approximations and lets the
// network compute F'(x) and F''(x) analytically via the chain rule.
#ifndef PP_ACTIVATION_H
#define PP_ACTIVATION_H

#include <functional>
#include <string>

namespace pp {

// Bundles a scalar activation sigma and its derivatives/antiderivative.
// All five fields must be filled; the ANN assumes they are consistent.
struct Activation {
    std::string name;
    std::function<double(double)> f;          // sigma(x)
    std::function<double(double)> df;         // sigma'(x)
    std::function<double(double)> ddf;        // sigma''(x)
    std::function<double(double)> dddf;       // sigma'''(x) (needed for F'' chain-rule terms)
    std::function<double(double)> antideriv;  // integral of sigma from 0 to x
};

// Returns the Gaussian wavelet activation: sigma(x) = x * exp(-x^2).
// Derivative: sigma'(x) = (1 - 2x^2) * exp(-x^2).
Activation gaussian_wavelet_activation();

// Returns the Gaussian activation: sigma(x) = exp(-x^2).
// Derivative: sigma'(x) = -2x * exp(-x^2).
Activation gaussian_activation();

}  // namespace pp

#endif  // PP_ACTIVATION_H
