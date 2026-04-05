// "problems.h" header file.
// ODE right-hand sides, problem solvers, and data writer declarations.
#pragma once

#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "vector.h"

namespace pp {

// Right-hand side for u'' = -u written as first-order system.
vector harmonic_rhs(double x, vector y);

// Right-hand side for Lotka-Volterra predator-prey model.
vector lotka_volterra_rhs(double x, vector y);

// Right-hand side generator for relativistic orbit equation.
std::function<vector(double, vector)> relativistic_rhs(double eps);

// Right-hand side for planar three-body dynamics.
vector three_body_rhs(double t, vector z);

// Solves the harmonic oscillator IVP.
std::tuple<std::vector<double>, std::vector<vector>> solve_oscillator();

// Solves the Lotka-Volterra IVP.
std::tuple<std::vector<double>, std::vector<vector>> solve_lotka_volterra();

// Solves the relativistic orbit IVP for given initial conditions.
std::tuple<std::vector<double>, std::vector<vector>> solve_relativistic_orbit(double eps, double u0, double up0);

// Solves the figure-8 three-body IVP.
std::tuple<std::vector<double>, std::vector<vector>> solve_three_body();

// Writes generic x/y data where y has two components.
void write_pair_data(const std::string& filename,
        const std::vector<double>& x, const std::vector<vector>& y);

// Writes orbit data as (phi, u, x, y).
void write_orbit_data(const std::string& filename,
        const std::vector<double>& phi, const std::vector<vector>& y);

// Writes figure-8 data as (t, x1, y1, x2, y2, x3, y3).
void write_three_body_data(const std::string& filename,
        const std::vector<double>& t, const std::vector<vector>& z);

} // namespace pp
