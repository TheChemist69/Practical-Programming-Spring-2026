// "rkstep.cc" implementation file.
// Embedded Runge-Kutta ODE integrator and adaptive-step-size driver.

#include "rkstep.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace pp {

// Advances the solution of dy/dx = f(x,y) by step h
// Uses an embedded Runge-Kutta 4(5) rule (Dormand-Prince) to estimate both y(x+h) and the error.
std::tuple<vector, vector> rkstep45(std::function<vector(double, vector)> F,
	double x, vector y, double h) {
	vector k1 = F(x, y);
	vector k2 = F(x + h / 5.0, y + k1 * (h * 1.0 / 5.0));
	vector k3 = F(x + 3.0 * h / 10.0,
		y + k1 * (h * 3.0 / 40.0) + k2 * (h * 9.0 / 40.0));
	vector k4 = F(x + 3.0 * h / 5.0,
		y + k1 * (h * 3.0 / 10.0) + k2 * (h * -9.0 / 10.0) + k3 * (h * 6.0 / 5.0));
	vector k5 = F(x + h,
		y + k1 * (h * -11.0 / 54.0) + k2 * (h * 5.0 / 2.0)
		  + k3 * (h * -70.0 / 27.0) + k4 * (h * 35.0 / 27.0));
	vector k6 = F(x + 7.0 * h / 8.0,
		y + k1 * (h * 1631.0 / 55296.0) + k2 * (h * 175.0 / 512.0)
		  + k3 * (h * 575.0 / 13824.0) + k4 * (h * 44275.0 / 110592.0)
		  + k5 * (h * 253.0 / 4096.0));

	vector yh = y + k1 * (h * 37.0 / 378.0) + k3 * (h * 250.0 / 621.0)
			+ k4 * (h * 125.0 / 594.0) + k6 * (h * 512.0 / 1771.0);
	vector y4 = y + k1 * (h * 2825.0 / 27648.0) + k3 * (h * 18575.0 / 48384.0)
			+ k4 * (h * 13525.0 / 55296.0) + k5 * (h * 277.0 / 14336.0)
			+ k6 * (h * 1.0 / 4.0);

	return {yh, yh - y4};
}

// Adaptive step-size driver routines advances solution from initial point a to final point b
// Automatically scales h up and down to restrict error relative to eps and acc bounds.
std::tuple<std::vector<double>, std::vector<vector>> driver(
	std::function<vector(double, vector)> F,
	double a, double b,
	vector yinit,
	double h,
	double acc,
	double eps,
	double hmax) {
	double x = a;
	vector y = yinit;
	std::vector<double> xlist{ x };
	std::vector<vector> ylist{ y };

	// Ensure the initial step points toward the target endpoint.
	if ((b > a && h < 0.0) || (b < a && h > 0.0)) h = -h;

	// Optional hard cap on step-size for better sampling of smooth trajectories.
	if (std::isfinite(hmax) && hmax > 0.0 && std::abs(h) > hmax) {
		h = std::copysign(hmax, h);
	}

	while (true) {
		if ((b >= a && x >= b) || (b < a && x <= b)) return {xlist, ylist};
		if ((b >= a && x + h > b) || (b < a && x + h < b)) h = b - x;

		auto [yh, dy] = rkstep45(F, x, y, h);
		double tol = (acc + eps * yh.norm()) * std::sqrt(std::abs(h) / std::abs(b - a));
		double err = dy.norm();

		if (err <= tol) {
			x += h;
			y = yh;
			xlist.push_back(x);
			ylist.push_back(y);
		}

		if (err > 0.0) {
			double factor = std::pow(tol / err, 0.2) * 0.95;
			h *= std::min(factor, 2.0);
		} else {
			h *= 2.0;
		}

		if (std::isfinite(hmax) && hmax > 0.0) {
			h = std::copysign(std::min(std::abs(h), hmax), h);
		}

		if (std::abs(h) < 10.0 * std::numeric_limits<double>::epsilon()) {
			std::fprintf(stderr, "adaptive driver stalled near x = %.16g\n", x);
			return {xlist, ylist};
		}
	}
}

} // namespace pp
