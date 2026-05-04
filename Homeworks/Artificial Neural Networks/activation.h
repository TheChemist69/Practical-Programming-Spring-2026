// Activation function bundle for the ANN.
#ifndef PP_ACTIVATION_H
#define PP_ACTIVATION_H

#include <functional>
#include <string>

namespace pp {

// Holds the activation function and its derivatives/antiderivative.
struct Activation {
	std::string name;
	std::function<double(double)> f;
	std::function<double(double)> df;
	std::function<double(double)> ddf;
	std::function<double(double)> dddf;
	std::function<double(double)> antideriv;
};

// Gaussian wavelet activation: f(x) = x * exp(-x^2).
Activation gaussian_wavelet_activation();

// Gaussian activation: f(x) = exp(-x^2).
Activation gaussian_activation();

}  // namespace pp

#endif  // PP_ACTIVATION_H
