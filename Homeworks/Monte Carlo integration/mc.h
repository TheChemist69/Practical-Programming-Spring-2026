// "mc.h" header file.
// Plain and stratified Monte Carlo integration declarations.
#pragma once

#include <cstddef>
#include <functional>
#include <random>

#include "vector.h"

namespace pp {

enum class MCStatus {
    success,
    invalid_input,
    non_finite_evaluation,
    max_depth_reached
};

struct MCResult {
    double value = 0.0;
    double estimated_error = 0.0;
    std::size_t evaluations = 0;
    MCStatus status = MCStatus::success;
};

struct StratifiedOptions {
    int nmin = 64;
    int max_depth = 64;
};

MCResult plain_mc(
        const std::function<double(const vector&)>& f,
        const vector& a,
        const vector& b,
        int n,
        std::mt19937_64& rng);

MCResult stratified_mc(
        const std::function<double(const vector&)>& f,
        const vector& a,
        const vector& b,
        int n,
        const StratifiedOptions& options,
        std::mt19937_64& rng);

const char* mc_status_cstr(MCStatus status);

} // namespace pp
