// "main.cc" -- Exam project driver.
// Adaptive integration of complex-valued functions along straight-line contours.
// This file ONLY orchestrates: it calls the routines in cintegrator / bessel,
// prints the report as a three-part exercise, runs an automated self-test, and
// writes the .data files consumed by the gnuplot scripts.
#include <array>
#include <cmath>
#include <complex>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <numbers>
#include <utility>
#include <vector>

#include "bessel.h"
#include "cintegrator.h"

namespace pp {
namespace {

using std::numbers::pi;
using namespace std::complex_literals;  // the imaginary unit 1i

// All console output is routed through emit(), which mirrors every line to the
// open log file as well.  Running the program therefore always produces Out.txt
// containing the complete printed report, no shell redirection needed.
// (A type-safe std::format wrapper would need a variadic template; this project
// deliberately uses no templates, so the C-style variadic tee is kept.)
FILE* g_logfile = nullptr;

void emit(const char* fmt, ...) {
    std::va_list ap;
    va_start(ap, fmt);
    std::vprintf(fmt, ap);  // to the console
    va_end(ap);
    if (g_logfile != nullptr) {
        va_start(ap, fmt);
        std::vfprintf(g_logfile, fmt, ap);  // and to Out.txt
        va_end(ap);
    }
}

// RAII owner of the log file: opens it, points emit() at it, closes on scope exit.
struct LogScope {
    explicit LogScope(const char* path) { g_logfile = std::fopen(path, "w"); }
    ~LogScope() {
        if (g_logfile != nullptr) {
            std::fclose(g_logfile);
            g_logfile = nullptr;
        }
    }
    LogScope(const LogScope&) = delete;
    LogScope& operator=(const LogScope&) = delete;
};

// ---- presentation helpers ---------------------------------------------------

const char* kRule =
    "──────────────────────────────────────────────────────────────";

// The title block printed once at the top of the report.
void title_block() {
    emit("\n  Adaptive Integration of Complex-Valued Functions\n");
    emit("  Exam project · Practical Programming 2026 · Nikolaos\n");
    emit("  %s\n", kRule);
}

// A major part header, e.g. part("Part A", "...", "6 pts").
void part(const char* label, const char* title, const char* pts) {
    emit("\n\n%s\n", kRule);
    emit("  %s · %s   [%s]\n", label, title, pts);
    emit("%s\n", kRule);
}

// A sub-topic heading inside a part (e.g. the B1/B2/B3 groups).
void topic(const char* s) { emit("\n  %s\n", s); }

// A single demonstration heading.
void demo(const char* s) { emit("\n  • %s\n", s); }

// ---- shared integrands (used by both the demos and the self-test) -----------

// Vertices of a regular polygon inscribed in the unit circle.
std::vector<cplx> unit_circle(int sides) {
    std::vector<cplx> v;
    v.reserve(sides);
    for (int k = 0; k < sides; k++) {
        const double t = 2.0 * pi * k / sides;
        v.push_back(cplx(std::cos(t), std::sin(t)));
    }
    return v;
}

// Contour integrand for ∫₀^{2π} dθ/(a + b cosθ) on |z| = 1, via z = e^{iθ}:
//   dθ = dz/(i z),  cosθ = (z + 1/z)/2.
std::function<cplx(cplx)> trig_integrand(double a, double b) {
    return [a, b](cplx z) -> cplx {
        return 1.0 / ((a + b * (z + 1.0 / z) * 0.5) * 1i * z);
    };
}

// f(z) = 1/(1 + z²), used for the rectangular-contour real-line integral.
cplx one_over_1_plus_z2(cplx z) { return 1.0 / (1.0 + z * z); }

// =============================================================================
// PART A: the complex-valued adaptive contour integrator
// =============================================================================
void run_part_a() {
    part("Part A", "A complex-valued adaptive contour integrator", "6 pts");

    demo("Analytic antiderivatives:  ∫ f(z) dz = F(b) − F(a)");
    {
        const cplx a(0.0, 0.0), b(1.0, 1.0);

        // integral of e^z from a to b = e^b - e^a
        const auto fexp = [](cplx z) { return std::exp(z); };
        const CResult rexp = integrate_segment(fexp, a, b);
        const cplx exp_exact = std::exp(b) - std::exp(a);
        emit("      ∫ e^z dz   over  0 → 1+i\n");
        emit("          computed  (% .10f, % .10f)\n", rexp.value.real(), rexp.value.imag());
        emit("          exact     (% .10f, % .10f)\n", exp_exact.real(), exp_exact.imag());
        emit("          error     %.1e   (%zu evaluations)\n",
             std::abs(rexp.value - exp_exact), rexp.evaluations);

        // integral of z^3 from a to b = (b^4 - a^4)/4
        const auto fcube = [](cplx z) { return z * z * z; };
        const CResult rcube = integrate_segment(fcube, a, b);
        const cplx cube_exact = (std::pow(b, 4) - std::pow(a, 4)) / 4.0;
        emit("      ∫ z^3 dz   over  0 → 1+i\n");
        emit("          computed  (% .10f, % .10f)\n", rcube.value.real(), rcube.value.imag());
        emit("          exact     (% .10f, % .10f)\n", cube_exact.real(), cube_exact.imag());
        emit("          error     %.1e\n", std::abs(rcube.value - cube_exact));
    }

    demo("Path independence for an analytic integrand");
    {
        const cplx a(-1.0, 0.0), b(1.0, 0.0);
        const auto f = [](cplx z) { return std::exp(z); };

        // Direct segment a -> b.
        const CResult direct = integrate_segment(f, a, b);
        // Detour a -> i -> b (two segments through the upper half-plane).
        const cplx via(0.0, 1.0);
        const CResult detour = integrate_segment(f, a, via);
        const CResult leg2 = integrate_segment(f, via, b);
        const cplx detour_val = detour.value + leg2.value;

        emit("      direct   −1 → 1        (% .10f, % .10f)\n",
             direct.value.real(), direct.value.imag());
        emit("      detour   −1 → i → 1    (% .10f, % .10f)\n",
             detour_val.real(), detour_val.imag());
        emit("      difference             %.1e\n", std::abs(direct.value - detour_val));

        std::ofstream out("partA_path_independence.data");
        out << "# label real imag\n";
        out << "direct " << direct.value.real() << " " << direct.value.imag() << "\n";
        out << "detour " << detour_val.real() << " " << detour_val.imag() << "\n";
    }

    demo("Residue theorem on a square contour around the origin");
    {
        // f(z) = 1/z has residue 1 at 0  ->  oint dz/z = 2 pi i.
        const auto finv = [](cplx z) { return 1.0 / z; };
        const std::vector<cplx> square = {cplx(1, -1), cplx(1, 1), cplx(-1, 1), cplx(-1, -1)};
        const CResult r = integrate_polygon(finv, square);
        emit("      ∮ dz/z       (% .10f, % .10f)\n", r.value.real(), r.value.imag());
        emit("      2πi exact    (% .10f, % .10f)   error %.1e\n",
             0.0, 2.0 * pi, std::abs(r.value - 2.0 * pi * 1i));

        // f(z) = z^2 is analytic  ->  closed integral = 0.
        const auto fsq = [](cplx z) { return z * z; };
        const CResult r0 = integrate_polygon(fsq, square);
        emit("      ∮ z^2 dz     |value| %.1e   (analytic ⇒ 0)\n", std::abs(r0.value));
    }

    demo("Work vs accuracy: cost of the adaptive rule as the tolerance tightens");
    {
        // ∫ e^z dz over 0 → 1+i has the exact value e^(1+i) − 1, so the achieved
        // error is known exactly.  Sweeping the tolerance shows the work–precision
        // trade-off: evaluations grow as the requested accuracy tightens.
        const cplx a(0.0, 0.0), b(1.0, 1.0);
        const auto fexp = [](cplx z) { return std::exp(z); };
        const cplx exact = std::exp(b) - std::exp(a);

        emit("          acc        evaluations     achieved error\n");
        std::ofstream out("work_precision.data");
        out << "# evaluations error acc\n";
        const double eps_machine = 2.220446049250313e-16;  // double precision
        for (double acc : {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14}) {
            COptions opt;
            opt.acc = acc;
            opt.eps = acc;
            const CResult r = integrate_segment(fexp, a, b, opt);
            const double err = std::abs(r.value - exact);
            emit("        %.0e        %7zu        %.2e\n", acc, r.evaluations, err);
            // For the log-scale plot, floor the error at machine epsilon so the
            // round-off plateau (where the true error is effectively 0) stays
            // visible: those points sit on the floor while evaluations keep rising.
            out << r.evaluations << " " << std::max(err, eps_machine) << " " << acc << "\n";
        }
    }
}

// =============================================================================
// PART B: real integrals evaluated by contour methods
// =============================================================================
void run_part_b() {
    part("Part B", "Real integrals evaluated by contour methods", "3 pts");

    // ---- B1: Bessel functions via contour integration ----------------------
    topic("B1 · Bessel functions Jₙ(x) by contour integration");

    demo("Contour integral vs power-series reference  (x = 2.5)");
    {
        const double x = 2.5;
        emit("          n         contour            reference          error\n");
        for (int n = 0; n <= 3; n++) {
            const double jc = bessel_contour(n, x, 6, 1.0);
            const double jr = bessel_reference(n, x);
            emit("          %d    % .12f    % .12f    %.1e\n", n, jc, jr, std::abs(jc - jr));
        }
    }

    demo("Independence of the polygon shape  (J₀ at x = 2.5)");
    {
        const double x = 2.5;
        const double jref = bessel_reference(0, x);
        std::ofstream conv("bessel_polygon_convergence.data");
        conv << "# sides error\n";
        for (int sides : {3, 4, 6, 8, 12, 24}) {
            const double jc = bessel_contour(0, x, sides, 1.0);
            emit("          %2d-gon   % .12f    error %.1e\n", sides, jc, std::abs(jc - jref));
            conv << sides << " " << std::abs(jc - jref) << "\n";
        }
    }

    demo("Recurrence  Jₙ₋₁(x) + Jₙ₊₁(x) = (2n/x) Jₙ(x)");
    {
        const double x = 4.0;
        const int n = 2;
        const double lhs = bessel_contour(n - 1, x, 6) + bessel_contour(n + 1, x, 6);
        const double rhs = (2.0 * n / x) * bessel_contour(n, x, 6);
        emit("      x = 4,  n = 2:   LHS % .10f    RHS % .10f    diff %.1e\n",
             lhs, rhs, std::abs(lhs - rhs));
    }

    // Data file: J_0,1,2 over a range of x for plotting.
    {
        std::ofstream out("bessel_contour.data");
        out << "# x J0_contour J0_ref J1_contour J1_ref J2_contour J2_ref\n";
        for (int i = 0; i <= 200; i++) {
            const double x = 0.0 + 12.0 * i / 200.0;
            out << x;
            for (int n = 0; n <= 2; n++) {
                out << " " << bessel_contour(n, x, 6, 1.0) << " " << bessel_reference(n, x);
            }
            out << "\n";
        }
    }

    // ---- B2: trigonometric integrals over the unit circle ------------------
    topic("B2 · Trigonometric integrals via the unit-circle substitution");
    const std::vector<cplx> circle = unit_circle(128);

    demo("∫₀^{2π} dθ/(a + b cosθ) = 2π/√(a²−b²),  via a contour on |z| = 1");
    {
        emit("           a     b        contour            exact            error\n");
        const std::array<std::pair<double, double>, 3> cases = {{{2.0, 1.0}, {5.0, 4.0}, {3.0, 2.0}}};
        for (const auto& [a, b] : cases) {
            const CResult r = integrate_polygon(trig_integrand(a, b), circle);
            const double exact = 2.0 * pi / std::sqrt(a * a - b * b);
            emit("         %4.1f  %4.1f    % .10f    % .10f    %.1e\n",
                 a, b, r.value.real(), exact, std::abs(r.value.real() - exact));
        }

        // Figure data: fix a = 2 and sweep b in [0, 1.9].
        std::ofstream out("real_trig.data");
        out << "# b numeric exact\n";
        for (int i = 0; i <= 95; i++) {
            const double a = 2.0;
            const double b = 1.9 * i / 95.0;
            const CResult r = integrate_polygon(trig_integrand(a, b), circle);
            const double exact = 2.0 * pi / std::sqrt(a * a - b * b);
            out << b << " " << r.value.real() << " " << exact << "\n";
        }
    }

    demo("Cost as b → a, the pole approaching the contour  (a = 2)");
    {
        emit("            b        contour          error      evaluations\n");
        for (const double b : {1.5, 1.9, 1.95, 1.99}) {
            const double a = 2.0;
            const CResult r = integrate_polygon(trig_integrand(a, b), circle);
            const double exact = 2.0 * pi / std::sqrt(a * a - b * b);
            emit("         %6.3f    % .8f    %.1e    %zu\n",
                 b, r.value.real(), std::abs(r.value.real() - exact), r.evaluations);
        }
    }

    // ---- B3: a real-line integral via a rectangular contour ----------------
    topic("B3 · A real-line integral via a closed rectangular contour");

    // f(z) = 1/(1+z^2) has a pole at z = i with residue 1/(2i).  A rectangle
    // [-R,R] x [0,H] taken counter-clockwise encloses it, so the closed integral
    // equals 2πi·Res = π exactly.  The real integral ∫_{-R}^{R} dx/(1+x^2) is the
    // bottom edge; it tends to π because the other three edges vanish as R,H → ∞.

    demo("∮ around a rectangle enclosing z = i equals 2πi·Res = π");
    {
        const double R = 20.0, H = 20.0;
        const std::vector<cplx> rect = {cplx(-R, 0), cplx(R, 0), cplx(R, H), cplx(-R, H)};
        const CResult loop = integrate_polygon(one_over_1_plus_z2, rect);
        emit("      ∮ rectangle (R=H=%.0f)   (% .10f, % .10f)\n",
             R, loop.value.real(), loop.value.imag());
        emit("      2πi·Res(z=i) = π        % .10f    error %.1e\n",
             pi, std::abs(loop.value - pi));
    }

    demo("∫_{-∞}^{∞} dx/(1+x²) = π:  the bottom edge → π as the other edges vanish");
    {
        emit("            R      bottom edge ∫_{-R}^{R}     |π − bottom|\n");
        std::ofstream out("real_rectangle.data");
        out << "# R bottom edge_error\n";
        for (const double R : {2.0, 5.0, 10.0, 20.0, 40.0, 80.0}) {
            const CResult bottom = integrate_segment(one_over_1_plus_z2, cplx(-R, 0), cplx(R, 0));
            const double b = bottom.value.real();
            const double err = std::abs(pi - b);
            emit("          %4.0f      % .10f            %.2e\n", R, b, err);
            out << R << " " << b << " " << err << "\n";
        }
    }
}

// =============================================================================
// PART C: spectral accuracy on a periodic contour
// =============================================================================
void run_part_c() {
    part("Part C", "Spectral trapezoid vs a low-order adaptive rule", "1 pt");

    const int n = 3;
    const double x = 5.0;
    const double ref = bessel_reference(n, x);

    demo("Periodic trapezoid vs the (low-order) adaptive open-4 rule, both for J₃(5)");
    emit("      target   J₃(5) = %.12f\n", ref);

    // Spectral trapezoid: error vs number of samples (= evaluations).
    emit("\n      periodic trapezoid (evaluations = N)\n");
    {
        std::ofstream out("spectral_trapezoid.data");
        out << "# evaluations error\n";
        for (int N : {4, 6, 8, 10, 12, 16, 20, 24, 32}) {
            const double val = bessel_trapezoid(n, x, N);
            const double err = std::abs(val - ref);
            emit("          N = %2d      error %.2e\n", N, err);
            out << N << " " << err << "\n";
        }
    }

    // Adaptive open-4 on the same periodic integrand over [0, 2 pi].
    emit("\n      adaptive open-4 (varying tolerance)\n");
    {
        // A real integrand expressed as a complex one with zero imaginary part.
        const std::function<cplx(double)> g = [n, x](double th) -> cplx {
            return std::cos(n * th - x * std::sin(th));
        };
        std::ofstream out("spectral_adaptive.data");
        out << "# evaluations error\n";
        for (double acc : {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12}) {
            COptions opt;
            opt.acc = acc;
            opt.eps = acc;
            const CResult r = adaptive_unit(g, 0.0, 2.0 * pi, opt);
            const double val = r.value.real() / (2.0 * pi);
            const double err = std::abs(val - ref);
            emit("          evals %7zu    error %.2e\n", r.evaluations, err);
            out << r.evaluations << " " << err << "\n";
        }
    }
}

// =============================================================================
// Automated self-test: re-validate the headline results to fixed tolerances.
// =============================================================================
int g_checks = 0;
int g_failures = 0;

void check(const char* name, double got, double want, double tol) {
    ++g_checks;
    const double err = std::abs(got - want);
    const bool ok = (err <= tol);
    if (!ok) ++g_failures;
    emit("      [%s]  %-32s  err %.1e  (tol %.0e)\n", ok ? "PASS" : "FAIL", name, err, tol);
}

void run_checks() {
    part("Checks", "automated verification of the headline results", "self-test");
    emit("\n");

    // Part A: residue theorem.
    {
        const auto finv = [](cplx z) { return 1.0 / z; };
        const std::vector<cplx> square = {cplx(1, -1), cplx(1, 1), cplx(-1, 1), cplx(-1, -1)};
        const CResult r = integrate_polygon(finv, square);
        check("residue  Im(oint dz/z) = 2pi", r.value.imag(), 2.0 * pi, 1e-9);
        check("residue  Re(oint dz/z) = 0", r.value.real(), 0.0, 1e-9);
    }

    // Part B1: Bessel value and recurrence.
    check("Bessel   J0(2.5) vs series", bessel_contour(0, 2.5, 6), bessel_reference(0, 2.5), 1e-12);
    {
        const double x = 4.0;
        const double lhs = bessel_contour(1, x, 6) + bessel_contour(3, x, 6);
        const double rhs = (2.0 * 2 / x) * bessel_contour(2, x, 6);
        check("Bessel   recurrence (n=2, x=4)", lhs, rhs, 1e-9);
    }

    // Part B2: unit-circle trig integral.
    {
        const CResult r = integrate_polygon(trig_integrand(2.0, 1.0), unit_circle(128));
        check("trig     a=2,b=1 vs 2pi/sqrt3", r.value.real(), 2.0 * pi / std::sqrt(3.0), 1e-9);
    }

    // Part B3: rectangular-contour residue.
    {
        const double R = 20.0;
        const std::vector<cplx> rect = {cplx(-R, 0), cplx(R, 0), cplx(R, R), cplx(-R, R)};
        const CResult loop = integrate_polygon(one_over_1_plus_z2, rect);
        check("rect     oint dz/(1+z^2) = pi", loop.value.real(), pi, 1e-9);
    }

    // Part C: spectral trapezoid.
    check("trap     J3(5) at N=32", bessel_trapezoid(3, 5.0, 32), bessel_reference(3, 5.0), 1e-12);

    emit("\n      %d checks, %d passed, %d failed\n", g_checks, g_checks - g_failures, g_failures);
}

int run() {
    const LogScope log("Out.txt");  // RAII: mirrors all output to Out.txt

    title_block();
    run_part_a();
    run_part_b();
    run_part_c();
    run_checks();

    emit("\n%s\n", kRule);
    emit("  Full console report saved to Out.txt.  Self-test: %d/%d checks passed.\n",
         g_checks - g_failures, g_checks);
    emit("%s\n", kRule);

    return (g_failures == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

}  // namespace
}  // namespace pp

int main() { return pp::run(); }
