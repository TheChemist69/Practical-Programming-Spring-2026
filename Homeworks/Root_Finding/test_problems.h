// "test_problems.h" header file.
// Analytic test systems used in Part A and Part C.
#pragma once

#include "newton.h"

namespace pp {

VectorFunction simple_scalar_equation();
VectorFunction simple_two_dimensional_system();
VectorFunction rosenbrock_gradient();
VectorFunction himmelblau_gradient();

double rosenbrock_value(const vector& x);
double himmelblau_value(const vector& x);

} // namespace pp
