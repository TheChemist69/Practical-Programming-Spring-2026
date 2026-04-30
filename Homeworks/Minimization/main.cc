// main.cc
// Runs minimization tests: Rosenbrock, Himmelblau, and Higgs fit
// Part C: Compares forward and central difference schemes

#include "minimization.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace pp;

static double rosenbrock(const std::vector<double>& x) {
    double X = x[0];
    double Y = x[1];
    return (1.0 - X)*(1.0 - X) + 100.0*(Y - X*X)*(Y - X*X);
}

static double himmelblau(const std::vector<double>& x) {
    double X = x[0];
    double Y = x[1];
    return (X*X + Y - 11.0)*(X*X + Y - 11.0) + (X + Y*Y - 7.0)*(X + Y*Y - 7.0);
}

static double breit_wigner_func(double E, double m, double Gamma, double A) {
    double den = (E - m)*(E - m) + 0.25*Gamma*Gamma;
    return A / den;
}

int main(int argc, char** argv) {
    std::cout << std::fixed << std::setprecision(6);
    
    // Configure optimizer
    MinimizerOptions opts;
    opts.acc = 1e-8;
    opts.max_iters = 1000;
    opts.use_central = false;   // Forward differences
    opts.levenberg = 0.1;        // Moderate damping
    
    // Rosenbrock: f(x,y) = (1-x)^2 + 100(y-x^2)^2
    std::vector<double> x0 = {-1.2, 1.0};
    auto res_r = newton_minimize(rosenbrock, x0, opts);
    std::cout << "Rosenbrock: iterations=" << res_r.iterations << " converged=" << res_r.converged
              << " x=[" << res_r.x[0] << ", " << res_r.x[1] << "] f=" << res_r.value << "\n";

    // Himmelblau: f(x,y) = (x^2+y-11)^2 + (x+y^2-7)^2
    // Note: Has 4 local minima at (3,2), (-2.8,3.1), (-3.8,-3.3), (3.6,-1.8)
    x0 = {0.0, 0.0};
    MinimizerOptions opts_h = opts;
    opts_h.levenberg = 0.5;      // Stronger damping for saddle points
    auto res_h = newton_minimize(himmelblau, x0, opts_h);
    std::cout << "Himmelblau: iterations=" << res_h.iterations << " converged=" << res_h.converged
              << " x=[" << res_h.x[0] << ", " << res_h.x[1] << "] f=" << res_h.value << "\n";

    // ============================================
    // PART C: Comparison of forward vs central differences
    // ============================================
    std::cout << "\n============================================\n";
    std::cout << "PART C: Forward vs Central Differences\n";
    std::cout << "============================================\n\n";
    
    std::cout << "Testing on Rosenbrock from (-1.2, 1.0):\n\n";
    
    // Forward differences (O(h) error, N function evals per gradient)
    MinimizerOptions opts_fwd = opts;
    opts_fwd.use_central = false;
    x0 = {-1.2, 1.0};
    auto res_fwd = newton_minimize(rosenbrock, x0, opts_fwd);
    std::cout << "Forward differences (O(h)):\n";
    std::cout << "  Iterations: " << res_fwd.iterations << "\n";
    std::cout << "  Solution: x=[" << res_fwd.x[0] << ", " << res_fwd.x[1] << "]\n";
    std::cout << "  Value: f=" << res_fwd.value << "\n";
    std::cout << "  Converged: " << (res_fwd.converged ? "yes" : "no") << "\n";
    
    // Central differences (O(h²) error, 2N function evals per gradient)
    MinimizerOptions opts_cen = opts;
    opts_cen.use_central = true;
    x0 = {-1.2, 1.0};
    auto res_cen = newton_minimize(rosenbrock, x0, opts_cen);
    std::cout << "\nCentral differences (O(h²)):\n";
    std::cout << "  Iterations: " << res_cen.iterations << "\n";
    std::cout << "  Solution: x=[" << res_cen.x[0] << ", " << res_cen.x[1] << "]\n";
    std::cout << "  Value: f=" << res_cen.value << "\n";
    std::cout << "  Converged: " << (res_cen.converged ? "yes" : "no") << "\n";
    
    std::cout << "\nComparison Summary:\n";
    std::cout << "  Forward:  " << res_fwd.iterations << " iterations (faster, less accurate)\n";
    std::cout << "  Central:  " << res_cen.iterations << " iterations (slower, more accurate)\n";
    std::cout << "  Central requires 2x more function evaluations per gradient computation\n";
    std::cout << "  but provides O(h²) truncation error vs O(h) for forward differences\n";

    // Higgs fit: read data from stdin: E sigma error
    std::vector<double> E, sigma, err;
    double a, b, c;
    while (std::cin >> a >> b >> c) {
        E.push_back(a);
        sigma.push_back(b);
        err.push_back(c);
    }
    if (E.empty()) {
        std::cerr << "No Higgs data provided on stdin; skipping fit.\n";
        return 0;
    }

    auto D = [&](const std::vector<double>& p) -> double {
        double m = p[0], G = p[1], A = p[2];
        double s = 0.0;
        for (size_t i = 0; i < E.size(); ++i) {
            double fit = breit_wigner_func(E[i], m, G, A);
            double d = (fit - sigma[i]) / err[i];
            s += d * d;
        }
        return s;
    };

    std::vector<double> p0 = {125.0, 4.0, 200.0};
    MinimizerOptions opts3 = opts;
    opts3.acc = 1e-6;
    auto res_p = newton_minimize(D, p0, opts3);
    std::cout << "\nHiggs fit: iterations=" << res_p.iterations << " converged=" << res_p.converged << "\n";
    std::cout << " m = " << res_p.x[0] << "  Gamma = " << res_p.x[1] << "  A = " << res_p.x[2]
              << "  D=" << res_p.value << "\n";

    // write fit curve to fit.out and also copy data to data.out
    FILE* fout = fopen("fit.out", "w");
    if (fout) {
        double Emin = *std::min_element(E.begin(), E.end());
        double Emax = *std::max_element(E.begin(), E.end());
        for (double ee = Emin; ee <= Emax; ee += 0.5) {
            double fv = breit_wigner_func(ee, res_p.x[0], res_p.x[1], res_p.x[2]);
            fprintf(fout, "%g %g\n", ee, fv);
        }
        fclose(fout);
    }
    FILE* dout = fopen("data.out", "w");
    if (dout) {
        for (size_t i = 0; i < E.size(); ++i)
            fprintf(dout, "%g %g %g\n", E[i], sigma[i], err[i]);
        fclose(dout);
    }

    std::cout << "Wrote fit.out and data.out for plotting.\n";
    return 0;
}
