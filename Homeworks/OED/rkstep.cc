// "rkstep.cc" implementation file.
// Embedded Runge-Kutta ODE integrator (Dormand-Prince RK45) and
// adaptive-step-size driver.
//
// The Dormand-Prince method uses 6 function evaluations per step and
// produces two solutions of different orders:
//   yh (5th order) -- used as the accepted solution
//   y4 (4th order) -- used only for error estimation: err = ||yh - y4||
// This is the "FSAL" (First Same As Last) trick: k6 is reused as k1
// of the next step, saving one evaluation.  (Not exploited here for
// simplicity; the driver always calls F fresh.)
//
// Butcher tableau (Dormand-Prince):
//   c = [0, 1/5, 3/10, 3/5, 1, 7/8]
//   5th-order weights b:  [37/378, 0, 250/621, 125/594, 0, 512/1771]
//   4th-order weights b*: [2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 1/4]
#include "rkstep.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace pp {

// Takes one Dormand-Prince RK45 step of size h from (x, y).
// Returns {yh, error_estimate} where yh is the 5th-order solution and
// error_estimate = yh - y4 is the difference between the 5th- and 4th-order solutions.
std::tuple<vector, vector> rkstep45(std::function<vector(double, vector)> F,
	double x, vector y, double h) {
	// Stage evaluations following the Dormand-Prince Butcher tableau.
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

	// 5th-order solution (accepted value).
	vector yh = y + k1 * (h * 37.0 / 378.0) + k3 * (h * 250.0 / 621.0)
			+ k4 * (h * 125.0 / 594.0) + k6 * (h * 512.0 / 1771.0);
	// 4th-order solution (used only to form the error estimate).
	vector y4 = y + k1 * (h * 2825.0 / 27648.0) + k3 * (h * 18575.0 / 48384.0)
			+ k4 * (h * 13525.0 / 55296.0) + k5 * (h * 277.0 / 14336.0)
			+ k6 * (h * 1.0 / 4.0);

	return {yh, yh - y4};
}

// Adaptive step-size driver: integrates F from x=a to x=b.
// At each step, tol = (acc + eps*||yh||) * sqrt(|h|/|b-a|) is the local
// tolerance budget.  If err <= tol the step is accepted; in either case
// the step size is rescaled by (tol/err)^0.2 * 0.95 (safety factor).
// The 0.2 exponent is optimal for 5th-order methods: it corresponds to
// h_new ~ h_old * (tol/err)^{1/(p+1)} with p=4 (the order of the error estimate).
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
