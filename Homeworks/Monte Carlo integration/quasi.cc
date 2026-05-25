// "quasi.cc" implementation file.
// Low-discrepancy Monte Carlo integration with two-sequence error estimate.
//
// Core idea:
// - Build deterministic low-discrepancy points from radical-inverse maps
//   (Halton-style sequence in multiple prime bases).
// - Run two sequence variants and compare them to estimate uncertainty.

#include "quasi.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace pp {

namespace {

bool valid_box(const vector& a, const vector& b) {
    // Same geometric validity criteria used by pseudo-random methods.
    if (a.size() <= 0 || a.size() != b.size()) {
        return false;
    }
    for (int i = 0; i < a.size(); i++) {
        if (!(b[i] > a[i])) {
            return false;
        }
    }
    return true;
}

double box_volume(const vector& a, const vector& b) {
    // Hyper-rectangle volume used to convert average value into integral.
    double volume = 1.0;
    for (int i = 0; i < a.size(); i++) {
        volume *= (b[i] - a[i]);
    }
    return volume;
}

double radical_inverse(std::uint64_t index, int base) {
    // Reflect the base-b digits of index across the decimal point.
    // Example in base 2: index 5 (101b) -> 0.101b = 0.625.
    double value = 0.0;
    double factor = 1.0 / static_cast<double>(base);

    while (index > 0) {
        value += factor * static_cast<double>(index % static_cast<std::uint64_t>(base));
        index /= static_cast<std::uint64_t>(base);
        factor /= static_cast<double>(base);
    }

    return value;
}

std::vector<int> first_primes(int n) {
    // Simple incremental prime generation; fast enough for small dimensions.
    std::vector<int> primes;
    if (n <= 0) {
        return primes;
    }

    primes.reserve(static_cast<std::size_t>(n));
    int candidate = 2;
    while (static_cast<int>(primes.size()) < n) {
        bool prime = true;
        for (int p : primes) {
            if (p * p > candidate) {
                break;
            }
            if (candidate % p == 0) {
                prime = false;
                break;
            }
        }
        if (prime) {
            primes.push_back(candidate);
        }
        candidate++;
    }

    return primes;
}

double halton_integral_estimate(
        const std::function<double(const vector&)>& f,
        const vector& a,
        const vector& b,
        int n,
        std::size_t start_index,
        const std::vector<int>& bases,
        MCStatus& status,
        std::size_t& evals) {
    // Build one deterministic estimate using a fixed set of Halton bases.
    const int dim = a.size();
    vector x(dim);

    double sum = 0.0;
    status = MCStatus::success;
    evals = 0;

    for (int i = 0; i < n; i++) {
        // Sequence index for this sample.
        const std::uint64_t idx = static_cast<std::uint64_t>(start_index + static_cast<std::size_t>(i));
        for (int k = 0; k < dim; k++) {
            // Convert [0,1) low-discrepancy coordinate to [a_k,b_k].
            const double u = radical_inverse(idx, bases[static_cast<std::size_t>(k)]);
            x[k] = a[k] + u * (b[k] - a[k]);
        }

        const double fx = f(x);
        evals++;

        if (!std::isfinite(fx)) {
            status = MCStatus::non_finite_evaluation;
            return std::numeric_limits<double>::quiet_NaN();
        }

        sum += fx;
    }

    return sum / n * box_volume(a, b);
}

} // namespace

MCResult quasi_mc_two_sequences(
        const std::function<double(const vector&)>& f,
        const vector& a,
        const vector& b,
        int n,
        const QuasiOptions& options) {
    // Reject malformed domain or impossible sequence parameters.
    if (!valid_box(a, b) || n <= 0 || options.start_index == 0 || options.first_base_offset < 0 || options.second_base_offset < 0) {
        return MCResult{
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::infinity(),
                0,
                MCStatus::invalid_input,
        };
    }

    const int dim = a.size();
    // Build enough primes for both base-offset windows.
    const int max_offset = (options.first_base_offset > options.second_base_offset)
            ? options.first_base_offset
            : options.second_base_offset;
    const std::vector<int> primes = first_primes(max_offset + dim + 1);

    std::vector<int> bases_a(static_cast<std::size_t>(dim));
    std::vector<int> bases_b(static_cast<std::size_t>(dim));
    for (int k = 0; k < dim; k++) {
        // Sequence A and B use different base blocks.
        bases_a[static_cast<std::size_t>(k)] = primes[static_cast<std::size_t>(options.first_base_offset + k)];
        bases_b[static_cast<std::size_t>(k)] = primes[static_cast<std::size_t>(options.second_base_offset + k)];
    }

    MCStatus status_a = MCStatus::success;
    MCStatus status_b = MCStatus::success;
    std::size_t evals_a = 0;
    std::size_t evals_b = 0;

    const double q1 = halton_integral_estimate(f, a, b, n, options.start_index, bases_a, status_a, evals_a);
    const double q2 = halton_integral_estimate(f, a, b, n, options.start_index, bases_b, status_b, evals_b);

    MCStatus status = MCStatus::success;
    if (status_a != MCStatus::success) {
        status = status_a;
    }
    if (status_b != MCStatus::success) {
        status = status_b;
    }

    if (status != MCStatus::success || !std::isfinite(q1) || !std::isfinite(q2)) {
        return MCResult{
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::infinity(),
                evals_a + evals_b,
                MCStatus::non_finite_evaluation,
        };
    }

    return MCResult{
            // Best estimate: average of both deterministic sequences.
            0.5 * (q1 + q2),
            // Error proxy: disagreement between the two sequence estimates.
            std::abs(q1 - q2),
            evals_a + evals_b,
            MCStatus::success,
    };
}

} // namespace pp
