// "newton.h" header file.
// Newton root solver with numerical Jacobian and configurable line search.
//
// Algorithm overview:
//   At each iteration, the Jacobian J = df/dx is built by forward finite
//   differences (one extra f-call per column).  The Newton step dx is found
//   by solving J*dx = -f(x) via QR decomposition.  A line search then
//   finds a step-size alpha in (0,1] such that the acceptance criterion is met.
//
//   Part A uses simple backtracking: accept as soon as ||f(x+alpha*dx)|| < ||f(x)||.
//   Part C uses quadratic interpolation for alpha, which requires fewer f-calls
//   and uses the stricter Armijo condition to measure improvement.
#pragma once

#include <functional>

#include "matrix.h"
#include "vector.h"

namespace pp {

// A vector-valued function R^n -> R^n.
using VectorFunction = std::function<vector(const vector&)>;

// Selects the line-search algorithm used to find the step size alpha.
enum class LineSearchType {
    backtracking,             // halve alpha until ||f(z)|| < ||f(x)|| (Part A)
    quadratic_interpolation   // fit a quadratic to phi(alpha) and jump to its minimum (Part C)
};

// Outcome of a newton_solve call.
enum class NewtonStatus {
    converged,          // ||f(x)|| < acc
    step_too_small,     // accepted step < finite-difference step scale (spec stopping criterion)
    max_iterations,     // ran out of iterations without converging
    singular_jacobian,  // QR solve failed (rank-deficient Jacobian)
    invalid_input       // bad arguments (empty x0, non-positive acc, etc.)
};

// User-tunable parameters for the Newton solver.
struct NewtonOptions {
    double acc = 1e-8;       // convergence threshold: stop when ||f(x)|| < acc
    double alpha_min = 1e-3; // minimum accepted step-size; the solver gives up shrinking below this
    int max_iter = 100;      // maximum number of Newton iterations

    // Optional user-supplied finite-difference step vector (length n).
    // If empty (the default), the solver computes safe default steps from x.
    vector dx;

    // Controls the default step size when dx is not supplied:
    //   true:  dx_i = max(|x_i|, 1) * 2^-26  (safe even when x_i = 0)
    //   false: dx_i = |x_i| * 2^-26           (proportional to x_i)
    bool default_dx_use_max_abs = true;

    LineSearchType line_search = LineSearchType::backtracking;
};

// Diagnostic counters collected during a solve.
struct NewtonStats {
    int iterations = 0;             // number of Newton steps taken
    int f_evaluations = 0;          // total calls to the user's f
    int jacobian_evaluations = 0;   // number of times the Jacobian was assembled
    int line_search_iterations = 0; // cumulative line-search refinement steps

    // The Jacobian matrix is allocated once and reused in-place each iteration
    // (not reallocated).  This counter should be 1 after a successful solve.
    int jacobian_allocations = 0;
};

// Output bundle returned by newton_solve.
struct NewtonResult {
    vector x;                                // best iterate found
    vector fx;                               // f(x) at the returned x
    NewtonStatus status = NewtonStatus::invalid_input;
    NewtonStats stats;
};

// Computes a safe finite-difference step vector for the given x.
// use_max_abs = true  => dx_i = max(|x_i|, 1) * sqrt(eps_machine)
// use_max_abs = false => dx_i = |x_i| * sqrt(eps_machine)
vector finite_difference_step(const vector& x, bool use_max_abs);

// Fills the Jacobian matrix J by forward finite differences.
// J_{ij} = (f_i(x + dx_j*e_j) - f_i(x)) / dx_j
// fx must be a pre-computed f(x) to avoid a redundant evaluation.
// If stats is non-null, increments f_evaluations and jacobian_evaluations.
void numerical_jacobian(
        const VectorFunction& f,
        const vector& x,
        const vector& fx,
        const vector& dx,
        matrix& J,
        NewtonStats* stats = nullptr);

// Finds a root of f starting from x0.
// Returns when ||f(x)|| < options.acc, when the step becomes smaller than
// the finite-difference scale, or when max_iter is reached.
NewtonResult newton_solve(const VectorFunction& f, const vector& x0, const NewtonOptions& options = {});

// Human-readable labels for status and line-search type (for output/logging).
const char* newton_status_cstr(NewtonStatus status);
const char* line_search_cstr(LineSearchType type);

} // namespace pp
