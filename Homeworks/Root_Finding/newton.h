// "newton.h" header file.
// Newton root solver with numerical Jacobian and line search.
#pragma once

#include <functional>

#include "matrix.h"
#include "vector.h"

namespace pp {

using VectorFunction = std::function<vector(const vector&)>;

enum class LineSearchType {
    backtracking,
    quadratic_interpolation
};

enum class NewtonStatus {
    converged,
    step_too_small,
    max_iterations,
    singular_jacobian,
    invalid_input
};

struct NewtonOptions {
    double acc = 1e-8;
    double alpha_min = 1e-3;
    int max_iter = 100;

    // Optional user-supplied finite-difference step vector.
    // If empty, the solver computes a default step from x.
    vector dx;

    // Default finite-difference scheme when dx is not supplied.
    // true:  dx_i = max(|x_i|,1)*2^-26
    // false: dx_i = |x_i|*2^-26
    bool default_dx_use_max_abs = true;

    LineSearchType line_search = LineSearchType::backtracking;
};

struct NewtonStats {
    int iterations = 0;
    int f_evaluations = 0;
    int jacobian_evaluations = 0;
    int line_search_iterations = 0;

    // Requirement-C evidence: Jacobian is allocated once per solve.
    int jacobian_allocations = 0;
};

struct NewtonResult {
    vector x;
    vector fx;
    NewtonStatus status = NewtonStatus::invalid_input;
    NewtonStats stats;
};

vector finite_difference_step(const vector& x, bool use_max_abs);
void numerical_jacobian(
        const VectorFunction& f,
        const vector& x,
        const vector& fx,
        const vector& dx,
        matrix& J,
        NewtonStats* stats = nullptr);

NewtonResult newton_solve(const VectorFunction& f, const vector& x0, const NewtonOptions& options = {});

const char* newton_status_cstr(NewtonStatus status);
const char* line_search_cstr(LineSearchType type);

} // namespace pp
