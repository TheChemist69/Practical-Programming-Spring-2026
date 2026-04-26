// "mc.cc" implementation file.
// Plain and recursive stratified Monte Carlo integration.
//
// Implementation overview:
// - plain_mc performs uniform random sampling in a box [a,b].
// - stratified_mc recursively partitions [a,b] and allocates samples
//   where pilot statistics indicate larger variance.

#include "mc.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace pp {

namespace {

// Compact running statistics used for pilot samples.
// We keep sum and sum of squares to estimate variance without
// storing all sampled values.
struct Stats {
    int count = 0;
    double sum = 0.0;
    double sum2 = 0.0;

    void add(double x) {
        count++;
        sum += x;
        sum2 += x * x;
    }
};

double variance_from_stats(const Stats& s) {
    // For empty stats we return zero variance so callers can keep logic simple.
    if (s.count <= 0) {
        return 0.0;
    }
    const double mean = s.sum / s.count;
    double var = s.sum2 / s.count - mean * mean;
    if (var < 0.0) {
        var = 0.0;
    }
    return var;
}

double box_volume(const vector& a, const vector& b) {
    // Hyper-rectangle volume = product of side lengths.
    double volume = 1.0;
    for (int i = 0; i < a.size(); i++) {
        volume *= (b[i] - a[i]);
    }
    return volume;
}

bool valid_box(const vector& a, const vector& b) {
    // Input is valid only when dimensions match and each upper bound is larger.
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

MCStatus merge_status(MCStatus a, MCStatus b) {
    // Return the most severe status observed in recursive subcalls.
    if (a == MCStatus::non_finite_evaluation || b == MCStatus::non_finite_evaluation) {
        return MCStatus::non_finite_evaluation;
    }
    if (a == MCStatus::invalid_input || b == MCStatus::invalid_input) {
        return MCStatus::invalid_input;
    }
    if (a == MCStatus::max_depth_reached || b == MCStatus::max_depth_reached) {
        return MCStatus::max_depth_reached;
    }
    return MCStatus::success;
}

MCResult stratified_recursive(
        const std::function<double(const vector&)>& f,
        const vector& a,
        const vector& b,
        int n,
        const StratifiedOptions& options,
        std::mt19937_64& rng,
        int depth) {
    // A branch may receive zero budget after variance-based splitting.
    if (n <= 0) {
        return MCResult{0.0, 0.0, 0, MCStatus::success};
    }

    // Very deep splitting can collapse an interval due to finite precision.
    // Treat the resulting measure-zero box as a harmless exhausted branch.
    for (int i = 0; i < a.size(); i++) {
        if (!(b[i] > a[i])) {
            return MCResult{0.0, 0.0, 0, MCStatus::max_depth_reached};
        }
    }

    if (depth >= options.max_depth) {
        // Stop refining and finish with plain sampling in the current sub-box.
        MCResult fallback = plain_mc(f, a, b, n, rng);
        fallback.status = MCStatus::max_depth_reached;
        return fallback;
    }

    if (n < options.nmin) {
        // Not enough budget to run a meaningful pilot phase.
        return plain_mc(f, a, b, n, rng);
    }

    const int dim = a.size();
    const int pilot = options.nmin;
    const int remaining = n - pilot;

    std::uniform_real_distribution<double> u01(0.0, 1.0);

    vector mid(dim);
    for (int i = 0; i < dim; i++) {
        // Candidate split location for each axis.
        mid[i] = 0.5 * (a[i] + b[i]);
    }

    Stats total;
    std::vector<Stats> left(dim);
    std::vector<Stats> right(dim);

    vector x(dim);
    for (int i = 0; i < pilot; i++) {
        // Draw one uniform point in the current box.
        for (int k = 0; k < dim; k++) {
            x[k] = a[k] + u01(rng) * (b[k] - a[k]);
        }

        const double fx = f(x);
        if (!std::isfinite(fx)) {
            return MCResult{
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::infinity(),
                    static_cast<std::size_t>(i + 1),
                    MCStatus::non_finite_evaluation,
            };
        }

        total.add(fx);
        for (int k = 0; k < dim; k++) {
            // For each dimension we record whether this sample lies
            // on the left or right side of the midpoint split.
            if (x[k] < mid[k]) {
                left[k].add(fx);
            } else {
                right[k].add(fx);
            }
        }
    }

    if (remaining <= 0) {
        // Pilot-only estimate if no samples remain for recursive refinement.
        const double mean = total.sum / pilot;
        const double sigma = std::sqrt(variance_from_stats(total));
        const double volume = box_volume(a, b);
        return MCResult{
                mean * volume,
                sigma * volume / std::sqrt(static_cast<double>(pilot)),
                static_cast<std::size_t>(pilot),
                MCStatus::success,
        };
    }

    int best_dim = 0;
    double best_measure = -1.0;
    for (int k = 0; k < dim; k++) {
        // Heuristic: split where left+right pilot variance is largest.
        const double measure = variance_from_stats(left[k]) + variance_from_stats(right[k]);
        if (measure > best_measure) {
            best_measure = measure;
            best_dim = k;
        }
    }

    double left_var = variance_from_stats(left[best_dim]);
    double right_var = variance_from_stats(right[best_dim]);

    if (left_var <= 0.0 && right_var <= 0.0) {
        // Flat pilot statistics: split budget evenly by using equal pseudo-variance.
        left_var = 1.0;
        right_var = 1.0;
    }

    int n_left = 0;
    int n_right = 0;

    if (remaining == 1) {
        // One sample left: send it to the noisier side.
        n_left = (left_var >= right_var) ? 1 : 0;
        n_right = 1 - n_left;
    } else {
        // Allocate in proportion to pilot variance.
        const double ratio = left_var / (left_var + right_var);
        n_left = static_cast<int>(std::lround(ratio * remaining));
        n_left = std::clamp(n_left, 1, remaining - 1);
        n_right = remaining - n_left;
    }

    vector left_a = a;
    vector left_b = b;
    vector right_a = a;
    vector right_b = b;

    left_b[best_dim] = mid[best_dim];
    right_a[best_dim] = mid[best_dim];

    // Recurse independently on both children and combine integrals.
    MCResult left_res = stratified_recursive(f, left_a, left_b, n_left, options, rng, depth + 1);
    MCResult right_res = stratified_recursive(f, right_a, right_b, n_right, options, rng, depth + 1);

    MCResult out;
    out.value = left_res.value + right_res.value;
    out.estimated_error = std::sqrt(
            left_res.estimated_error * left_res.estimated_error
            + right_res.estimated_error * right_res.estimated_error);
    out.evaluations = static_cast<std::size_t>(pilot) + left_res.evaluations + right_res.evaluations;
    out.status = merge_status(left_res.status, right_res.status);
    return out;
}

} // namespace

MCResult plain_mc(
        const std::function<double(const vector&)>& f,
        const vector& a,
        const vector& b,
        int n,
        std::mt19937_64& rng) {
    // Guard against malformed domain or invalid sample count.
    if (!valid_box(a, b) || n <= 0) {
        return MCResult{
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::infinity(),
                0,
                MCStatus::invalid_input,
        };
    }

    const int dim = a.size();
    const double volume = box_volume(a, b);

    std::uniform_real_distribution<double> u01(0.0, 1.0);
    vector x(dim);

    double sum = 0.0;
    double sum2 = 0.0;

    for (int i = 0; i < n; i++) {
        // Generate one random point uniformly in [a,b].
        for (int k = 0; k < dim; k++) {
            x[k] = a[k] + u01(rng) * (b[k] - a[k]);
        }

        const double fx = f(x);
        if (!std::isfinite(fx)) {
            return MCResult{
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::infinity(),
                    static_cast<std::size_t>(i + 1),
                    MCStatus::non_finite_evaluation,
            };
        }

        sum += fx;
        sum2 += fx * fx;
    }

    const double mean = sum / n;
    double sigma2 = sum2 / n - mean * mean;
    if (sigma2 < 0.0) {
        sigma2 = 0.0;
    }

    return MCResult{
            // Integral estimate = average f(x) times domain volume.
            mean * volume,
            // Standard error of the Monte Carlo mean times volume.
            std::sqrt(sigma2 / n) * volume,
            static_cast<std::size_t>(n),
            MCStatus::success,
    };
}

MCResult stratified_mc(
        const std::function<double(const vector&)>& f,
        const vector& a,
        const vector& b,
        int n,
        const StratifiedOptions& options,
        std::mt19937_64& rng) {
    // Validate user options before entering recursion.
    if (!valid_box(a, b) || n <= 0 || options.nmin <= 0 || options.max_depth <= 0) {
        return MCResult{
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::infinity(),
                0,
                MCStatus::invalid_input,
        };
    }

    return stratified_recursive(f, a, b, n, options, rng, 0);
}

const char* mc_status_cstr(MCStatus status) {
    // Small utility used in console logs and diagnostics.
    switch (status) {
        case MCStatus::success:
            return "success";
        case MCStatus::invalid_input:
            return "invalid_input";
        case MCStatus::non_finite_evaluation:
            return "non_finite_evaluation";
        case MCStatus::max_depth_reached:
            return "max_depth_reached";
    }
    return "unknown";
}

} // namespace pp
