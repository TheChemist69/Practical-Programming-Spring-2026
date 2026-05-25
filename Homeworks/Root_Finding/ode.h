// "ode.h" header file.
// Adaptive ODE integrator used by the hydrogen shooting method.
#pragma once

#include <cstddef>
#include <functional>
#include <vector>

#include "vector.h"

namespace pp {

using ODEFunction = std::function<vector(double, const vector&)>;

struct ODEOptions {
    double acc = 1e-6;
    double eps = 1e-6;
    double h_initial = 0.05;
    double h_min = 1e-8;
    int max_steps = 200000;
};

struct ODEStats {
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    std::size_t function_evaluations = 0;
    bool reached_end = true;
    bool step_underflow = false;
};

struct ODEResult {
    std::vector<double> xs;
    std::vector<vector> ys;
    vector y_end;
    ODEStats stats;
};

ODEResult integrate_ode_adaptive(
        const ODEFunction& f,
        double a,
        const vector& ya,
        double b,
        const ODEOptions& options,
        bool store_path);

} // namespace pp
