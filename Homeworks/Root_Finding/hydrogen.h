// "hydrogen.h" header file.
// Shooting-method utilities for the hydrogen radial equation.
#pragma once

#include <vector>

#include "ode.h"

namespace pp {

struct ShootingConfig {
    double rmin = 1e-3;
    double rmax = 8.0;

    ODEOptions ode_options;

    double e_min = -1.5;
    double e_max = -0.01;
    int scan_points = 400;

    double root_acc = 1e-8;
    int root_max_iter = 120;
};

struct ShootingResult {
    double energy = 0.0;
    double mismatch = 0.0;
    bool converged = false;
    int iterations = 0;
};

struct WavefunctionPoint {
    double r = 0.0;
    double f_numeric = 0.0;
    double f_exact = 0.0;
};

struct ConvergenceRow {
    double parameter = 0.0;
    double energy = 0.0;
    double abs_error = 0.0;
    bool converged = false;
};

double hydrogen_exact_energy_ground();
double hydrogen_exact_wavefunction(double r);

double hydrogen_boundary_mismatch(double energy, const ShootingConfig& cfg, ODEStats* ode_stats = nullptr);
ShootingResult find_ground_state_energy(const ShootingConfig& cfg);
std::vector<WavefunctionPoint> compute_wavefunction(const ShootingConfig& cfg, double energy);

} // namespace pp
