// "mc.h" header file.
// Public interface for pseudo-random Monte Carlo integration methods.
//
// This module exposes two estimators:
// 1) plain_mc: classic uniform random sampling in a hyper-rectangle.
// 2) stratified_mc: recursive variance-guided splitting of the domain.
//
// Both methods return a value estimate, an internal error estimate,
// the number of function evaluations, and a status code.
#pragma once

#include <cstddef>
#include <functional>
#include <random>

#include "vector.h"

namespace pp {

// Execution state of a Monte Carlo call.
// The status allows the caller to distinguish numerical failure
// from normal completion.
enum class MCStatus {
    success,
    invalid_input,
    non_finite_evaluation,
    max_depth_reached
};

// Unified result container used by all Monte Carlo methods.
// - value: estimated integral.
// - estimated_error: method-specific uncertainty estimate.
// - evaluations: number of integrand calls actually performed.
// - status: completion state.
struct MCResult {
    double value = 0.0;
    double estimated_error = 0.0;
    std::size_t evaluations = 0;
    MCStatus status = MCStatus::success;
};

// Controls for recursive stratified sampling.
// - nmin: number of pilot samples used before deciding a split.
// - max_depth: recursion guard against endless refinement.
struct StratifiedOptions {
    int nmin = 64;
    int max_depth = 64;
};

// Plain Monte Carlo on the axis-aligned box [a,b].
// The integrand is sampled uniformly, and the returned error is
// derived from the sample variance.
MCResult plain_mc(
        const std::function<double(const vector&)>& f,
        const vector& a,
        const vector& b,
        int n,
        std::mt19937_64& rng);

// Recursive stratified Monte Carlo on [a,b].
// The algorithm spends a pilot budget, estimates directional
// variability, then splits and allocates the remaining budget
// recursively where variance appears larger.
MCResult stratified_mc(
        const std::function<double(const vector&)>& f,
        const vector& a,
        const vector& b,
        int n,
        const StratifiedOptions& options,
        std::mt19937_64& rng);

// Convert a status enum to a human-readable string for logs.
const char* mc_status_cstr(MCStatus status);

} // namespace pp
