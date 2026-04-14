// "quasi.h" header file.
// Low-discrepancy Monte Carlo integration declarations.
#pragma once

#include <cstddef>
#include <functional>

#include "mc.h"
#include "vector.h"

namespace pp {

struct QuasiOptions {
    std::size_t start_index = 1;
    int first_base_offset = 0;
    int second_base_offset = 32;
};

MCResult quasi_mc_two_sequences(
        const std::function<double(const vector&)>& f,
        const vector& a,
        const vector& b,
        int n,
        const QuasiOptions& options = {});

} // namespace pp
