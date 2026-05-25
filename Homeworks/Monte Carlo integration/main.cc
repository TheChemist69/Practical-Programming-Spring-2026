// "main.cc" implementation file.
// Orchestration for Monte Carlo integration homework tasks.
//
// This file does not implement integration algorithms directly.
// Instead, it coordinates experiments, computes diagnostic metrics,
// prints readable summaries, and writes machine-readable .data files.

#include "mc.h"
#include "problems.h"
#include "quasi.h"
#include "reporting.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

namespace pp {

namespace {

// Fit y = c * N^p on log-log scale by linear regression.
// Input vectors are interpreted as:
// - ns[i]    -> N
// - errors[i]-> observed error at N
// Returned value is slope p.
double fit_loglog_slope(const std::vector<int>& ns, const std::vector<double>& errors) {
    double sx = 0.0;
    double sy = 0.0;
    double sxx = 0.0;
    double sxy = 0.0;
    int m = 0;

    for (std::size_t i = 0; i < ns.size() && i < errors.size(); i++) {
        if (ns[i] <= 0 || errors[i] <= 0.0) {
            // Log transform requires strictly positive values.
            continue;
        }
        const double x = std::log(static_cast<double>(ns[i]));
        const double y = std::log(errors[i]);
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
        m++;
    }

    if (m < 2) {
        // Not enough valid points to fit a line.
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double denom = m * sxx - sx * sx;
    if (std::abs(denom) < 1e-16) {
        // Degenerate x distribution (all points nearly identical).
        return std::numeric_limits<double>::quiet_NaN();
    }

    return (m * sxy - sx * sy) / denom;
}

void run_two_dimensional_examples() {
    // Small sanity check on two standard 2D problems.
    // This gives a quick confidence signal before larger scaling sweeps.
    std::puts("=== Two-dimensional test integrals ===");

    const MCProblem circle = unit_circle_area_problem();
    const MCProblem smooth = smooth_gaussian_2d_problem();

    const int n = 150000;
    // Fixed seed for reproducible homework output.
    std::mt19937_64 rng(2026041401ull);

    MCResult circle_res = plain_mc(circle.integrand, circle.a, circle.b, n, rng);
    MCResult smooth_res = plain_mc(smooth.integrand, smooth.a, smooth.b, n, rng);

    std::vector<TwoDExampleRow> rows;
    rows.push_back(TwoDExampleRow{
            "plain",
            circle.name,
            n,
            circle_res.value,
            circle.exact,
            std::abs(circle_res.value - circle.exact),
            circle_res.estimated_error,
    });
    rows.push_back(TwoDExampleRow{
            "plain",
            smooth.name,
            n,
            smooth_res.value,
            smooth.exact,
            std::abs(smooth_res.value - smooth.exact),
            smooth_res.estimated_error,
    });

    write_two_dim_examples("two_dim_examples.data", rows);

    std::printf("  %-20s N=%d  I=% .12g exact=% .12g |err|=%.3e est=%.3e\n",
            circle.name.c_str(),
            n,
            circle_res.value,
            circle.exact,
            std::abs(circle_res.value - circle.exact),
            circle_res.estimated_error);

    std::printf("  %-20s N=%d  I=% .12g exact=% .12g |err|=%.3e est=%.3e\n\n",
            smooth.name.c_str(),
            n,
            smooth_res.value,
            smooth.exact,
            std::abs(smooth_res.value - smooth.exact),
            smooth_res.estimated_error);
}

void run_plain_scaling() {
    // Study how plain Monte Carlo error decreases with N.
    std::puts("=== Plain Monte Carlo scaling study ===");

    const MCProblem problem = smooth_gaussian_2d_problem();
    const std::vector<int> ns{200, 500, 1000, 2000, 5000, 10000, 20000, 50000};

    std::vector<PlainScalingRow> rows;
    rows.reserve(ns.size());

    std::vector<double> actual_errors;
    actual_errors.reserve(ns.size());

    for (int n : ns) {
        // Per-N deterministic seed keeps runs reproducible while still
        // varying random streams between sample counts.
        std::mt19937_64 rng(910000ull + static_cast<unsigned long long>(n));
        MCResult res = plain_mc(problem.integrand, problem.a, problem.b, n, rng);
        const double actual = std::abs(res.value - problem.exact);

        rows.push_back(PlainScalingRow{
                n,
                res.estimated_error,
                actual,
                res.value,
                problem.exact,
            // Multiplying by sqrt(N) helps inspect the expected O(N^-1/2)
            // convergence behavior of plain Monte Carlo.
                res.estimated_error * std::sqrt(static_cast<double>(n)),
                actual * std::sqrt(static_cast<double>(n)),
        });
        actual_errors.push_back(actual);

        std::printf("  N=%-7d value=% .12g |err|=% .3e est=% .3e status=%s\n",
                n,
                res.value,
                actual,
                res.estimated_error,
                mc_status_cstr(res.status));
    }

    write_plain_scaling("plain_scaling.data", rows);

    const double slope = fit_loglog_slope(ns, actual_errors);
    std::printf("  log-log slope(actual_error vs N): %.4f (expected about -0.5)\n\n", slope);
}

void run_quasi_scaling() {
    // Compare pseudo-random MC against quasi-random MC on the same problem.
    std::puts("=== Quasi-random vs pseudo-random scaling ===");

    const MCProblem problem = smooth_gaussian_2d_problem();
    const std::vector<int> ns{200, 500, 1000, 2000, 5000, 10000, 20000, 50000};

    std::vector<QuasiScalingRow> rows;
    rows.reserve(ns.size());

    std::vector<double> plain_actual;
    std::vector<double> quasi_actual;
    plain_actual.reserve(ns.size());
    quasi_actual.reserve(ns.size());

    for (int n : ns) {
        // Separate pseudo-random seed sequence for this study.
        std::mt19937_64 rng(920000ull + static_cast<unsigned long long>(n));
        MCResult plain = plain_mc(problem.integrand, problem.a, problem.b, n, rng);
        MCResult quasi = quasi_mc_two_sequences(problem.integrand, problem.a, problem.b, n);

        const double plain_err = std::abs(plain.value - problem.exact);
        const double quasi_err = std::abs(quasi.value - problem.exact);

        rows.push_back(QuasiScalingRow{
                n,
                plain_err,
                quasi_err,
                plain.estimated_error,
                quasi.estimated_error,
                plain.value,
                quasi.value,
                problem.exact,
        });

        plain_actual.push_back(plain_err);
        quasi_actual.push_back(quasi_err);

        std::printf("  N=%-7d plain|err|=% .3e quasi|err|=% .3e\n", n, plain_err, quasi_err);
    }

    write_quasi_scaling("quasi_scaling.data", rows);

    const double slope_plain = fit_loglog_slope(ns, plain_actual);
    const double slope_quasi = fit_loglog_slope(ns, quasi_actual);
    std::printf("  slope plain = %.4f, slope quasi = %.4f\n\n", slope_plain, slope_quasi);
}

void run_stratified_scaling() {
    // Compare recursive stratified sampling against plain sampling.
    std::puts("=== Stratified vs plain scaling ===");

    const MCProblem problem = smooth_gaussian_2d_problem();
    const std::vector<int> ns{200, 500, 1000, 2000, 5000, 10000, 20000, 50000};

    StratifiedOptions strat_opts;
    // Baseline recursive controls used in this project.
    strat_opts.nmin = 64;
    strat_opts.max_depth = 64;

    std::vector<StratifiedScalingRow> rows;
    rows.reserve(ns.size());

    for (int n : ns) {
        // Independent seeds avoid accidental cross-correlation.
        std::mt19937_64 plain_rng(930000ull + static_cast<unsigned long long>(n));
        std::mt19937_64 strat_rng(940000ull + static_cast<unsigned long long>(n));

        MCResult plain = plain_mc(problem.integrand, problem.a, problem.b, n, plain_rng);
        MCResult strat = stratified_mc(problem.integrand, problem.a, problem.b, n, strat_opts, strat_rng);

        const double plain_err = std::abs(plain.value - problem.exact);
        const double strat_err = std::abs(strat.value - problem.exact);

        rows.push_back(StratifiedScalingRow{
                n,
                plain_err,
                strat_err,
                plain.estimated_error,
                strat.estimated_error,
                static_cast<double>(plain.evaluations),
                static_cast<double>(strat.evaluations),
                problem.exact,
        });

        std::printf("  N=%-7d plain|err|=% .3e strat|err|=% .3e calls_plain=%zu calls_strat=%zu\n",
                n,
                plain_err,
                strat_err,
                plain.evaluations,
                strat.evaluations);
    }

    write_stratified_scaling("stratified_scaling.data", rows);
    std::puts("");
}

void run_difficult_integral() {
    // Hard benchmark with near-singular behavior.
    // We report convergence trends but avoid over-interpreting asymptotic slope.
    std::puts("=== Difficult singular integral (not used for scaling fits) ===");

    const MCProblem hard = difficult_singular_3d_problem();
    const std::vector<int> ns{1000, 3000, 10000, 30000};

    StratifiedOptions strat_opts;
    // Slightly larger pilot and depth budget to better resolve local features.
    strat_opts.nmin = 96;
    strat_opts.max_depth = 72;

    std::vector<DifficultIntegralRow> rows;
    rows.reserve(ns.size());

    for (int n : ns) {
        std::mt19937_64 plain_rng(950000ull + static_cast<unsigned long long>(n));
        std::mt19937_64 strat_rng(960000ull + static_cast<unsigned long long>(n));

        MCResult plain = plain_mc(hard.integrand, hard.a, hard.b, n, plain_rng);
        MCResult quasi = quasi_mc_two_sequences(hard.integrand, hard.a, hard.b, n);
        MCResult strat = stratified_mc(hard.integrand, hard.a, hard.b, n, strat_opts, strat_rng);

        rows.push_back(DifficultIntegralRow{
                n,
                plain.value,
                quasi.value,
                strat.value,
                std::abs(plain.value - hard.exact),
                std::abs(quasi.value - hard.exact),
                std::abs(strat.value - hard.exact),
                hard.exact,
        });

        std::printf("  N=%-6d plain=% .9f quasi=% .9f strat=% .9f\n",
                n,
                plain.value,
                quasi.value,
                strat.value);
    }

    write_difficult_integral("difficult_integral.data", rows);

    const DifficultIntegralRow& last = rows.back();
    std::printf("  exact = %.16f\n", hard.exact);
    std::printf("  largest-N errors: plain=%.3e quasi=%.3e strat=%.3e\n\n",
            last.actual_error_plain,
            last.actual_error_quasi,
            last.actual_error_stratified);
}

} // namespace

int run() {
    // Execute all experiment sections in a stable order so output files
    // and console summaries are consistent across runs.
    run_lcg_comparison();          // Part A: LCG vs Mersenne Twister
    run_two_dimensional_examples();
    run_plain_scaling();
    run_quasi_scaling();
    run_stratified_scaling();
    run_difficult_integral();

    std::puts("Data files written:");
    std::puts("  lcg_comparison.data");
    std::puts("  two_dim_examples.data");
    std::puts("  plain_scaling.data");
    std::puts("  quasi_scaling.data");
    std::puts("  stratified_scaling.data");
    std::puts("  difficult_integral.data");

    return EXIT_SUCCESS;
}

// Part A: comparison of the LCG against std::mt19937_64.
// The spec asks students to compare three generators; here we run both on the
// same smooth 2D problem and tabulate errors at increasing N.
void run_lcg_comparison() {
    std::puts("=== Part A: LCG vs Mersenne Twister comparison ===");

    const MCProblem problem = smooth_gaussian_2d_problem();
    const std::vector<int> ns{200, 500, 1000, 2000, 5000, 10000, 20000, 50000};

    std::printf("  %-7s  %-12s  %-12s\n", "N", "|err| LCG", "|err| MT64");

    FILE* fp = std::fopen("lcg_comparison.data", "w");
    if (fp) std::fprintf(fp, "# N err_lcg err_mt64\n");

    for (int n : ns) {
        // Different seeds but same N so the comparison is fair.
        Lcg lcg(static_cast<uint32_t>(n) ^ 0xDEADBEEFu);
        std::mt19937_64 mt(970000ull + static_cast<unsigned long long>(n));

        MCResult lcg_res = plain_mc_lcg(problem.integrand, problem.a, problem.b, n, lcg);
        MCResult mt_res  = plain_mc(problem.integrand, problem.a, problem.b, n, mt);

        const double lcg_err = std::abs(lcg_res.value - problem.exact);
        const double mt_err  = std::abs(mt_res.value  - problem.exact);

        std::printf("  %-7d  %-12.3e  %-12.3e\n", n, lcg_err, mt_err);
        if (fp) std::fprintf(fp, "%d %.6e %.6e\n", n, lcg_err, mt_err);
    }

    if (fp) std::fclose(fp);
    std::puts("");
}

} // namespace pp

int main() {
    // Thin program entry point; real logic lives in pp::run().
    return pp::run();
}
