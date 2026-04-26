// "quasi.h" header file.
// Interface for low-discrepancy (Halton-based) Monte Carlo integration.
//
// The implementation runs two independent Halton-style sequences and
// uses their difference as a pragmatic uncertainty indicator.
#pragma once

#include <cstddef>
#include <functional>

#include "mc.h"
#include "vector.h"

namespace pp {

// Configuration for the two-sequence quasi-random estimator.
// - start_index: first Halton index (must be >= 1).
// - first_base_offset / second_base_offset: where to start choosing
//   prime bases for sequence A and sequence B.
// Different base offsets make the two sequences decorrelated enough
// to produce a useful error proxy.
struct QuasiOptions {
    std::size_t start_index = 1;
    int first_base_offset = 0;
    int second_base_offset = 32;
};

// Quasi-random Monte Carlo estimate on [a,b].
// Returns the mean of two sequence estimates and uses |q1-q2| as
// an empirical error estimate.
MCResult quasi_mc_two_sequences(
        const std::function<double(const vector&)>& f,
        const vector& a,
        const vector& b,
        int n,
        const QuasiOptions& options = {});

} // namespace pp
