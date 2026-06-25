// "rkstep.h" header file.
// Embedded Runge-Kutta stepper and adaptive driver declarations.
#pragma once

#include <functional>
#include <limits>
#include <tuple>
#include <vector>

#include "vector.h"

namespace pp {

// Advances y by one RK45 step and returns {y(x+h), estimated error}.
std::tuple<vector, vector> rkstep45(std::function<vector(double, vector)> F,
        double x, vector y, double h);

// Integrates from a to b using adaptive step-size control.
std::tuple<std::vector<double>, std::vector<vector>> driver(
        std::function<vector(double, vector)> F,
        double a, double b,
        vector yinit,
        double h,
        double acc,
        double eps,
        double hmax = std::numeric_limits<double>::infinity());

} // namespace pp
