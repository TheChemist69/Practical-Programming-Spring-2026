// "newton.cc" implementation file.
// Newton solver with numerical Jacobian and configurable line search.

#include "newton.h"

#include <algorithm>
#include <cmath>

#include "qr.h"

namespace pp {

namespace {

bool line_search_condition(double fz_norm, double fx_norm, double alpha) {
    return fz_norm < (1.0 - alpha * 0.5) * fx_norm;
}

double quadratic_alpha_update(double alpha, double fx_norm, double fz_norm) {
    // Quadratic model with phi(0)=||f(x)|| and approximate phi'(0)=-||f(x)||.
    // This gives a trial minimizer used as an improved alternative to alpha/=2.
    const double denom = 2.0 * (fz_norm - fx_norm + alpha * fx_norm);
    if (!(denom > 0.0)) {
        return 0.5 * alpha;
    }

    double trial = fx_norm * alpha * alpha / denom;

    // Keep progress conservative and monotone decreasing.
    const double lo = 0.1 * alpha;
    const double hi = 0.5 * alpha;
    trial = std::clamp(trial, lo, hi);

    if (!std::isfinite(trial) || trial <= 0.0) {
        trial = 0.5 * alpha;
    }
    return trial;
}

} // namespace

vector finite_difference_step(const vector& x, bool use_max_abs) {
    const double sqrt_machine_eps = std::ldexp(1.0, -26);
    vector dx(x.size());
    for (int i = 0; i < x.size(); i++) {
        const double xi_abs = std::abs(x[i]);
        const double scale = use_max_abs ? std::max(xi_abs, 1.0) : xi_abs;
        dx[i] = scale * sqrt_machine_eps;
        if (dx[i] == 0.0) {
            dx[i] = sqrt_machine_eps;
        }
    }
    return dx;
}

void numerical_jacobian(
        const VectorFunction& f,
        const vector& x,
        const vector& fx,
        const vector& dx,
        matrix& J,
        NewtonStats* stats) {
    const int n = x.size();

    // Requirement-C optimization: J is reused and overwritten in-place.
    for (int j = 0; j < n; j++) {
        vector x_shifted = x;
        x_shifted[j] += dx[j];

        const vector df = f(x_shifted) - fx;
        if (stats) {
            stats->f_evaluations++;
        }

        for (int i = 0; i < n; i++) {
            J(i, j) = df[i] / dx[j];
        }
    }

    if (stats) {
        stats->jacobian_evaluations++;
    }
}

NewtonResult newton_solve(const VectorFunction& f, const vector& x0, const NewtonOptions& options) {
    NewtonResult out;

    if (x0.size() <= 0 || options.max_iter <= 0 || options.acc <= 0.0 || options.alpha_min <= 0.0) {
        out.status = NewtonStatus::invalid_input;
        out.x = x0;
        out.fx = vector(x0.size());
        return out;
    }

    const int n = x0.size();
    out.x = x0;
    out.fx = f(out.x);
    out.stats.f_evaluations = 1;

    if (out.fx.size() != n) {
        out.status = NewtonStatus::invalid_input;
        return out;
    }

    matrix J(n, n);
    out.stats.jacobian_allocations = 1;

    for (int iter = 0; iter < options.max_iter; iter++) {
        const double fx_norm = out.fx.norm();
        if (fx_norm < options.acc) {
            // "iter" full Newton updates were already completed.
            out.stats.iterations = iter;
            out.status = NewtonStatus::converged;
            return out;
        }

        const vector dx_default = finite_difference_step(out.x, options.default_dx_use_max_abs);

        vector dx_fd;
        if (options.dx.size() == n) {
            dx_fd = options.dx;

            // Harden user-supplied finite-difference steps.
            // Non-finite or non-positive entries are replaced by safe defaults.
            for (int i = 0; i < n; i++) {
                if (!std::isfinite(dx_fd[i]) || dx_fd[i] <= 0.0) {
                    dx_fd[i] = dx_default[i];
                }
            }
        } else {
            dx_fd = dx_default;
        }

        numerical_jacobian(f, out.x, out.fx, dx_fd, J, &out.stats);

        bool solve_ok = true;
        const vector step = solve_linear_system_qr(J, -out.fx, &solve_ok);
        if (!solve_ok) {
            out.stats.iterations = iter;
            out.status = NewtonStatus::singular_jacobian;
            return out;
        }

        double alpha = 1.0;
        vector z = out.x;
        vector fz = out.fx;

        while (true) {
            z = out.x + alpha * step;
            fz = f(z);
            out.stats.f_evaluations++;

            const double fz_norm = fz.norm();
            if (line_search_condition(fz_norm, fx_norm, alpha)) {
                break;
            }

            if (alpha < options.alpha_min) {
                // Accept the poor step if alpha has already become too small.
                break;
            }

            if (options.line_search == LineSearchType::quadratic_interpolation) {
                alpha = quadratic_alpha_update(alpha, fx_norm, fz_norm);
            } else {
                alpha *= 0.5;
            }
            out.stats.line_search_iterations++;
        }

        const vector accepted_step = alpha * step;
        out.x = z;
        out.fx = fz;
        out.stats.iterations = iter + 1;

        // Assignment stop condition: step smaller than Jacobian finite-difference scale.
        if (accepted_step.norm() < dx_fd.norm()) {
            out.status = NewtonStatus::step_too_small;
            if (out.fx.norm() < options.acc) {
                out.status = NewtonStatus::converged;
            }
            return out;
        }
    }

    out.status = NewtonStatus::max_iterations;
    out.stats.iterations = options.max_iter;
    if (out.fx.norm() < options.acc) {
        out.status = NewtonStatus::converged;
    }
    return out;
}

const char* newton_status_cstr(NewtonStatus status) {
    switch (status) {
        case NewtonStatus::converged:
            return "converged";
        case NewtonStatus::step_too_small:
            return "step_too_small";
        case NewtonStatus::max_iterations:
            return "max_iterations";
        case NewtonStatus::singular_jacobian:
            return "singular_jacobian";
        case NewtonStatus::invalid_input:
            return "invalid_input";
    }
    return "unknown";
}

const char* line_search_cstr(LineSearchType type) {
    switch (type) {
        case LineSearchType::backtracking:
            return "backtracking";
        case LineSearchType::quadratic_interpolation:
            return "quadratic_interpolation";
    }
    return "unknown";
}

} // namespace pp
