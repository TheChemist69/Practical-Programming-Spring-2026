// Nikolaos Antonios Iakynthos Nemet
// This project was created using a combination of Claude Code and GitHub Copilot.

// "main.cc" implementation file.
// Entry orchestration for the ODE homework report runs.

#include "problems.h"

#include <cmath>
#include <cstdio>

namespace pp {

// ---------------------------------------------------------------
// Embedded Runge-Kutta stepper + adaptive driver
// ---------------------------------------------------------------
// This section demonstrates that the rkstep/driver implementation
// is working by solving representative IVPs requested in the text.
static void run_task_embedded_rk_driver() {
        std::puts("=== Task (6 points): Embedded RK + adaptive-step driver ===");
        std::puts("  Solving y' = f(x,y) with automatic step-size control.");

        // Debug/reference problem from the assignment: u'' = -u.
        // State form is y0' = y1, y1' = -y0.
        std::puts("\n  [A1] Harmonic oscillator test: u'' = -u");
        auto [x1, y1] = solve_oscillator();
        write_pair_data("oscillator.data", x1, y1);
        std::printf("    points: %zu, final |u-u_exact| = %.3e\n",
                x1.size(), std::abs(y1.back()[0] - std::cos(x1.back())));
        std::puts("    output: oscillator.data");

        // Also requested by the assignment as an interesting system:
        // Lotka-Volterra predator-prey model.
        std::puts("\n  [A2] Lotka-Volterra system (scipy-style benchmark)");
        auto [x2, y2] = solve_lotka_volterra();
        write_pair_data("lotka_volterra.data", x2, y2);
        std::printf("    points: %zu\n", x2.size());
        std::puts("    output: lotka_volterra.data");
}

// ------------------------------------------------------
// Relativistic orbital precession in GR
// ------------------------------------------------------
// Governing equation from assignment:
//   u''(phi) + u(phi) = 1 + eps * u(phi)^2
// Rewritten as first-order system in problems.cc.
static void run_task_relativistic_orbit() {
        std::puts("\n=== Task (3 points): Relativistic orbit precession ===");
        std::puts("  Equation: u'' + u = 1 + eps*u^2, with u = 1/r.");

        // Case B1: eps = 0, u(0)=1, u'(0)=0 -> Newtonian circular motion.
        std::puts("\n  [B1] Newtonian circular orbit: eps=0, u0=1, u0'=0");
        auto [phi1, orbit1] = solve_relativistic_orbit(0.0, 1.0, 0.0);
        write_orbit_data("orbit_newtonian_circle.data", phi1, orbit1);
        std::printf("    points: %zu\n", phi1.size());
        std::puts("    output: orbit_newtonian_circle.data");

        // Case B2: eps = 0, u(0)=1, u'(0)=-0.5 -> Newtonian ellipse.
        std::puts("\n  [B2] Newtonian elliptical orbit: eps=0, u0=1, u0'=-0.5");
        auto [phi2, orbit2] = solve_relativistic_orbit(0.0, 1.0, -0.5);
        write_orbit_data("orbit_newtonian_ellipse.data", phi2, orbit2);
        std::printf("    points: %zu\n", phi2.size());
        std::puts("    output: orbit_newtonian_ellipse.data");

        // Case B3: eps = 0.01, u(0)=1, u'(0)=-0.5 -> relativistic precession.
        std::puts("\n  [B3] Relativistic precession: eps=0.01, u0=1, u0'=-0.5");
        auto [phi3, orbit3] = solve_relativistic_orbit(0.01, 1.0, -0.5);
        write_orbit_data("orbit_relativistic.data", phi3, orbit3);
        std::printf("    points: %zu\n", phi3.size());
        std::puts("    output: orbit_relativistic.data");
}

// -------------------------------------------------------
// Figure-8 solution of the three-body ODE
// -------------------------------------------------------
// Integrates the planar 3-body system using the same adaptive
// ODE solver to reproduce the special periodic figure-8 motion.
static void run_task_three_body_figure8() {
        std::puts("\n=== Task (1 point): Three-body figure-8 solution ===");
        std::puts("  Integrating the 12-component planar 3-body state vector.");

        auto [t, z] = solve_three_body();
        write_three_body_data("figure8.data", t, z);
        std::printf("    points: %zu\n", t.size());
        std::puts("    output: figure8.data");
}

// Runs all tasks in assignment order and reports output files.
int run() {
        run_task_embedded_rk_driver();
        run_task_relativistic_orbit();
        run_task_three_body_figure8();

        std::puts("\nData files written:");
        std::puts("  oscillator.data");
        std::puts("  lotka_volterra.data");
        std::puts("  orbit_newtonian_circle.data");
        std::puts("  orbit_newtonian_ellipse.data");
        std::puts("  orbit_relativistic.data");
        std::puts("  figure8.data");
        return 0;
}

} // namespace pp

// Program entry point.
int main() {
        return pp::run();
}
