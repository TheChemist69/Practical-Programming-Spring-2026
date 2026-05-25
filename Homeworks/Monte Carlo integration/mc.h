// "mc.h" header file.
// Public interface for Monte Carlo integration methods.
//
// This module exposes:
// 1) Lcg: a quick-and-dirty Linear Congruential Generator (spec Part A).
//    Used to compare PRNG quality against the standard library Mersenne Twister.
// 2) plain_mc: classic uniform random sampling using std::mt19937_64.
// 3) plain_mc_lcg: same algorithm but driven by the LCG.
// 4) stratified_mc: recursive variance-guided splitting of the domain.
//
// All estimators return a value, an internal error estimate, evaluation count,
// and a status code.
#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <random>

#include "vector.h"

namespace pp {

// Quick-and-dirty Linear Congruential Generator (homework spec Part A).
// Parameters from Knuth: multiplier 1664525, increment 1013904223, modulus 2^32.
// Period = 2^32.  Quality is acceptable for low-dimensional MC but fails
// spectral tests and is unsuitable for cryptography or high-dimensional work.
// Used here to demonstrate convergence rate differences vs. std::mt19937_64.
class Lcg {
public:
    explicit Lcg(uint32_t seed = 1) : state_(seed) {}

    // Returns the next integer in [0, 2^32).
    uint32_t next() {
        state_ = 1664525u * state_ + 1013904223u;
        return state_;
    }

    // Returns a uniform sample in [0, 1).
    double uniform() {
        return static_cast<double>(next()) / 4294967296.0; // 2^32
    }

private:
    uint32_t state_;
};

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

// Plain Monte Carlo on [a,b] driven by std::mt19937_64 (high-quality PRNG).
// The integrand is sampled uniformly; the returned error is derived from
// the sample variance: sigma/sqrt(N) * volume.
MCResult plain_mc(
        const std::function<double(const vector&)>& f,
        const vector& a,
        const vector& b,
        int n,
        std::mt19937_64& rng);

// Same plain Monte Carlo algorithm but driven by the LCG (homework Part A).
// Useful for comparing convergence quality vs. the Mersenne Twister.
MCResult plain_mc_lcg(
        const std::function<double(const vector&)>& f,
        const vector& a,
        const vector& b,
        int n,
        Lcg& rng);

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
