// "hydrogen.h" header file.
// Shooting method for the hydrogen atom radial Schrödinger equation.
//
// The radial equation (atomic units, l=0 ground state) is:
//   f''(r) = -2(E + 1/r) f(r)
// This is cast as a first-order IVP from r = rmin with the asymptotic boundary
// condition f(rmin) ~ rmin - rmin^2, f'(rmin) ~ 1 - 2*rmin.
// The shooting method varies energy E until f(rmax) = 0 (mismatch M(E) = 0),
// giving the quantized ground-state energy E0 = -1/2 a.u.
#pragma once

#include <vector>

#include "ode.h"

namespace pp {

// Parameters controlling the shooting integration and root search.
struct ShootingConfig {
    double rmin = 1e-3;     // inner boundary (> 0; avoids the 1/r singularity)
    double rmax = 8.0;      // outer boundary (wavefunction must decay to ~0 by rmax)

    ODEOptions ode_options; // tolerances and step limits passed to the adaptive ODE driver

    double e_min = -1.5;    // lower bound of the energy scan range
    double e_max = -0.01;   // upper bound of the energy scan range
    int scan_points = 400;  // number of uniformly spaced energies probed when searching for a sign-change bracket

    double root_acc = 1e-8;  // bisection convergence criterion: stop when |M(E)| < root_acc
    int root_max_iter = 120; // maximum bisection iterations
};

// Outcome of a shooting root-finding call.
struct ShootingResult {
    double energy = 0.0;     // best energy estimate E0
    double mismatch = 0.0;   // residual M(energy) = f_E(rmax) at convergence
    bool converged = false;  // true if bisection terminated within root_acc
    int iterations = 0;      // number of bisection steps taken
};

// One sampled point along the radial wavefunction for comparison with the exact solution.
struct WavefunctionPoint {
    double r = 0.0;          // radial coordinate
    double f_numeric = 0.0;  // numerically integrated f(r)
    double f_exact = 0.0;    // analytic ground-state wavefunction: r*exp(-r)
};

// One row in a parameter-convergence table (vary rmin, rmax, or ODE tolerance).
struct ConvergenceRow {
    double parameter = 0.0;  // value of the varied parameter
    double energy = 0.0;     // computed ground-state energy for this parameter
    double abs_error = 0.0;  // |energy - exact| = |energy + 0.5|
    bool converged = false;  // whether the bisection converged
};

// Returns the analytically known hydrogen ground-state energy: E0 = -1/2 a.u.
double hydrogen_exact_energy_ground();

// Returns the analytic radial wavefunction for the 1s ground state: f(r) = r*exp(-r).
double hydrogen_exact_wavefunction(double r);

// Integrates the radial equation from rmin to rmax for the given energy and returns
// M(energy) = f_E(rmax).  A sign change in M(E) signals a bound-state energy between
// two scan points.  Optionally writes ODE statistics into ode_stats.
double hydrogen_boundary_mismatch(double energy, const ShootingConfig& cfg, ODEStats* ode_stats = nullptr);

// Finds the ground-state energy by scanning M(E) for a sign change and then
// bisecting to root_acc tolerance.
ShootingResult find_ground_state_energy(const ShootingConfig& cfg);

// Integrates the radial equation at the given energy and returns the trajectory
// downsampled to at most ~2000 points for plotting against the exact wavefunction.
std::vector<WavefunctionPoint> compute_wavefunction(const ShootingConfig& cfg, double energy);

} // namespace pp
