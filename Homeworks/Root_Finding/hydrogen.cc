// "hydrogen.cc" implementation file.
// Shooting method for hydrogen bound-state energy and wavefunction.

#include "hydrogen.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace pp {

namespace {

vector hydrogen_rhs(double energy, double r, const vector& y) {
    vector dydr(2);
    dydr[0] = y[1];
    dydr[1] = -2.0 * (energy + 1.0 / r) * y[0];
    return dydr;
}

vector hydrogen_initial_state(double rmin) {
    // Asymptotic boundary condition near r=0: f(r)=r-r^2.
    // Its derivative is f'(r)=1-2r.
    return vector{rmin - rmin * rmin, 1.0 - 2.0 * rmin};
}

} // namespace

double hydrogen_exact_energy_ground() {
    return -0.5;
}

double hydrogen_exact_wavefunction(double r) {
    return r * std::exp(-r);
}

double hydrogen_boundary_mismatch(double energy, const ShootingConfig& cfg, ODEStats* ode_stats) {
    if (!(cfg.rmin > 0.0) || !(cfg.rmax > cfg.rmin)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const vector y0 = hydrogen_initial_state(cfg.rmin);
    const ODEFunction f = [energy](double r, const vector& y) {
        return hydrogen_rhs(energy, r, y);
    };

    const ODEResult sol = integrate_ode_adaptive(f, cfg.rmin, y0, cfg.rmax, cfg.ode_options, false);

    if (ode_stats) {
        *ode_stats = sol.stats;
    }

    if (!sol.stats.reached_end || sol.y_end.size() != 2) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Shooting mismatch M(E) = F_E(rmax).
    return sol.y_end[0];
}

ShootingResult find_ground_state_energy(const ShootingConfig& cfg) {
    ShootingResult out;

    if (!(cfg.e_min < cfg.e_max) || cfg.scan_points < 2 || cfg.root_max_iter <= 0 || cfg.root_acc <= 0.0) {
        return out;
    }

    // First scan in energy to find the lowest sign-changing bracket.
    bool bracket_found = false;
    double a = 0.0;
    double b = 0.0;
    double fa = 0.0;

    double e_prev = cfg.e_min;
    double m_prev = hydrogen_boundary_mismatch(e_prev, cfg);

    for (int i = 1; i <= cfg.scan_points; i++) {
        const double t = static_cast<double>(i) / static_cast<double>(cfg.scan_points);
        const double e_cur = cfg.e_min + t * (cfg.e_max - cfg.e_min);
        const double m_cur = hydrogen_boundary_mismatch(e_cur, cfg);

        if (std::isfinite(m_prev) && std::isfinite(m_cur) && m_prev * m_cur <= 0.0) {
            a = e_prev;
            b = e_cur;
            fa = m_prev;
            bracket_found = true;
            break;
        }

        e_prev = e_cur;
        m_prev = m_cur;
    }

    if (!bracket_found) {
        return out;
    }

    // Bisection on M(E)=0 inside the bracket.
    double mid = 0.5 * (a + b);
    double fm = hydrogen_boundary_mismatch(mid, cfg);

    for (int iter = 0; iter < cfg.root_max_iter; iter++) {
        out.iterations = iter + 1;

        mid = 0.5 * (a + b);
        fm = hydrogen_boundary_mismatch(mid, cfg);

        if (!std::isfinite(fm)) {
            return out;
        }

        if (std::abs(fm) < cfg.root_acc || std::abs(b - a) < cfg.root_acc) {
            out.converged = true;
            out.energy = mid;
            out.mismatch = fm;
            return out;
        }

        if (fa * fm <= 0.0) {
            b = mid;
        } else {
            a = mid;
            fa = fm;
        }
    }

    out.energy = mid;
    out.mismatch = fm;
    return out;
}

std::vector<WavefunctionPoint> compute_wavefunction(const ShootingConfig& cfg, double energy) {
    const vector y0 = hydrogen_initial_state(cfg.rmin);
    const ODEFunction f = [energy](double r, const vector& y) {
        return hydrogen_rhs(energy, r, y);
    };

    const ODEResult sol = integrate_ode_adaptive(f, cfg.rmin, y0, cfg.rmax, cfg.ode_options, true);

    std::vector<WavefunctionPoint> rows;
    if (sol.xs.empty() || sol.ys.empty()) {
        return rows;
    }

    const std::size_t target_points = 2000;
    const std::size_t stride = std::max<std::size_t>(1, sol.xs.size() / target_points);
    rows.reserve(sol.xs.size() / stride + 2);

    for (std::size_t i = 0; i < sol.xs.size() && i < sol.ys.size(); i += stride) {
        rows.push_back(WavefunctionPoint{
                sol.xs[i],
                sol.ys[i][0],
                hydrogen_exact_wavefunction(sol.xs[i]),
        });
    }

    if (rows.back().r != sol.xs.back()) {
        rows.push_back(WavefunctionPoint{
                sol.xs.back(),
                sol.ys.back()[0],
                hydrogen_exact_wavefunction(sol.xs.back()),
        });
    }

    return rows;
}

} // namespace pp
