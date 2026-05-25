// "ode.cc" implementation file.
// Adaptive Runge-Kutta 1/2 integrator (Euler/midpoint pair).

#include "ode.h"

#include <algorithm>
#include <cmath>

namespace pp {

namespace {

void rkstep12(
        const ODEFunction& f,
        double x,
        const vector& y,
        double h,
        vector& y_next,
        vector& err,
        ODEStats& stats) {
    const vector k0 = f(x, y);
    stats.function_evaluations++;

    const vector y_mid = y + 0.5 * h * k0;
    const vector k1 = f(x + 0.5 * h, y_mid);
    stats.function_evaluations++;

    const vector y_euler = y + h * k0;
    y_next = y + h * k1;
    err = y_next - y_euler;
}

} // namespace

ODEResult integrate_ode_adaptive(
        const ODEFunction& f,
        double a,
        const vector& ya,
        double b,
        const ODEOptions& options,
        bool store_path) {
    ODEResult out;
    out.y_end = ya;

    if (ya.size() <= 0 || options.max_steps <= 0 || options.h_min <= 0.0) {
        out.stats.reached_end = false;
        return out;
    }

    if (a == b) {
        if (store_path) {
            out.xs.push_back(a);
            out.ys.push_back(ya);
        }
        return out;
    }

    const double direction = (b > a) ? 1.0 : -1.0;
    double h = std::abs(options.h_initial);
    if (h == 0.0) {
        h = std::abs(b - a) / 100.0;
    }
    h *= direction;

    double x = a;
    vector y = ya;

    if (store_path) {
        out.xs.push_back(x);
        out.ys.push_back(y);
    }

    for (int steps = 0; steps < options.max_steps; steps++) {
        if ((direction > 0.0 && x >= b) || (direction < 0.0 && x <= b)) {
            out.y_end = y;
            return out;
        }

        if ((direction > 0.0 && x + h > b) || (direction < 0.0 && x + h < b)) {
            h = b - x;
        }

        vector y_next(y.size());
        vector err(y.size());
        rkstep12(f, x, y, h, y_next, err, out.stats);

        const double err_norm = err.norm();
        const double y_norm = y_next.norm();

        // Local acceptance tolerance. A direct acc+eps*|y| criterion is
        // sufficient for this coursework and avoids over-restricting tiny steps.
        const double tol = options.acc + options.eps * y_norm;

        if (err_norm <= tol || err_norm == 0.0) {
            x += h;
            y = y_next;
            out.stats.accepted_steps++;

            if (store_path) {
                out.xs.push_back(x);
                out.ys.push_back(y);
            }
        } else {
            out.stats.rejected_steps++;
        }

        double factor = 2.0;
        if (err_norm > 0.0 && tol > 0.0) {
            factor = 0.95 * std::pow(tol / err_norm, 0.25);
            factor = std::clamp(factor, 0.2, 5.0);
        }

        h *= factor;

        if (std::abs(h) < options.h_min) {
            out.stats.step_underflow = true;
            out.stats.reached_end = false;
            out.y_end = y;
            return out;
        }
    }

    out.stats.reached_end = false;
    out.y_end = y;
    return out;
}

} // namespace pp
