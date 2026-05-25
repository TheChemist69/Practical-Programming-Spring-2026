// "main.cc" -- Ordinary Least-Squares Fit: ThX radioactive decay.
//
// Fits ln(activity) = ln(a) - lambda*t  to the 1902 Rutherford & Soddy
// measurements of ThX (224-Ra) decay.  Uses QR-decomposition-based
// ordinary least squares with proper uncertainty propagation.
#include <cstdio>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <limits>
#include "vector.h"
#include "matrix.h"
#include "lsfit.h"

int main() {
    std::printf("=== Ordinary Least-Squares Fit ===\n");
    std::printf("=== ThX (224-Ra) radioactive decay (Rutherford & Soddy, 1902) ===\n\n");

    // -- Experimental data ------------------------------------------
    // Time [days], activity [relative units], uncertainty in activity
    pp::vector t  = {1, 2, 3, 4, 6, 9, 10, 13, 15};
    pp::vector y  = {117, 100, 88, 72, 53, 29.5, 25.2, 15.2, 11.1};
    pp::vector dy = {6, 5, 4, 4, 4, 3, 3, 2, 2};

    int n = t.size();

    // -- Transform to log-space -------------------------------------
    // The decay law  y(t) = a*exp(-lambda*t)  becomes linear in logs:
    //   ln(y) = ln(a) - lambda*t = c_1*1 + c_2*t
    // Uncertainties propagate as  d(ln y) = dy / y  (first order).
    pp::vector lny(n);
    pp::vector dlny(n);
    for (int i = 0; i < n; i++) {
        lny[i]  = std::log(y[i]);
        dlny[i] = dy[i] / y[i];
    }

    // -- Basis functions for the log-linear model -------------------
    pp::BasisFunctions fs = {
        [](double /*t*/) { return 1.0; },   // f_1(t) = 1   -> coefficient c_1 = ln(a)
        [](double  t   ) { return t;   }    // f_2(t) = t   -> coefficient c_2 = -lambda
    };

    // -- Perform the fit --------------------------------------------
    pp::lsfit fit(fs, t, lny, dlny);

    // -- Extract physical parameters --------------------------------
    double ln_a      = fit.c[0];
    double neg_lambda = fit.c[1];            // c_2 = -lambda
    double a          = std::exp(ln_a);
    double lambda     = -neg_lambda;
    double half_life  = std::log(2.0) / lambda;

    // Uncertainties
    double d_ln_a      = fit.dc[0];
    double d_lambda    = fit.dc[1];          // |d(c_2)| = d_lambda
    // Error propagation:  T_half = ln2/lambda  =>  dT_half = (ln2/lambda^2)*d_lambda
    double d_half_life = std::log(2.0) / (lambda * lambda) * d_lambda;

    // -- Print results -----------------------------------------------
    std::printf("Fit: ln(y) = c1 + c2*t\n");
    std::printf("  c1 = ln(a) = %.6f +/- %.6f\n", ln_a, d_ln_a);
    std::printf("  c2 = -lam  = %.6f +/- %.6f\n", neg_lambda, d_lambda);
    std::printf("\n");
    std::printf("Derived quantities:\n");
    std::printf("  a      = %.4f\n", a);
    std::printf("  lambda = %.6f +/- %.6f  [1/day]\n", lambda, d_lambda);
    std::printf("  T_half = ln(2)/lambda = %.4f +/- %.4f  [days]\n",
                half_life, d_half_life);
    std::printf("  Modern value (224-Ra): T_half = 3.6319 days\n\n");

    std::printf("Covariance matrix of the fit coefficients:\n");
    fit.cov.print("  ");

    double rho = fit.cov[0, 1] / (fit.dc[0] * fit.dc[1]);
    std::printf("Correlation coefficient rho_12 = %.6f\n\n", rho);

    // -- Write data files for gnuplot --------------------------------
    constexpr double tmin = 0.0, tmax = 16.0;
    constexpr int    npts = 200;

    // 1. Original data with error bars
    {
        std::ofstream fout("data.txt");
        fout << "# t [days]   y [activity]   dy\n";
        for (int i = 0; i < n; i++)
            fout << t[i] << " " << y[i] << " " << dy[i] << "\n";
    }

    // 2. Best-fit curve and uncertainty bands (in original y-space)
    //    The uncertainty band is formed by evaluating the fit with all
    //    2^m combinations of c_k +/- dc_k and taking the envelope.
    {
        std::ofstream fout("fit.txt");
        fout << "# t   y_fit   y_upper   y_lower\n";
        for (int i = 0; i <= npts; i++) {
            double ti     = tmin + (tmax - tmin) * i / npts;
            double y_fit  = std::exp(fit(ti));

            // Enumerate +/- combinations for 2 coefficients (4 corners)
            double y_max = std::numeric_limits<double>::lowest();
            double y_min = std::numeric_limits<double>::max();
            for (int s1 = -1; s1 <= 1; s1 += 2) {
                for (int s2 = -1; s2 <= 1; s2 += 2) {
                    double ln_y_var = (fit.c[0] + s1 * fit.dc[0])
                                    + (fit.c[1] + s2 * fit.dc[1]) * ti;
                    double y_var = std::exp(ln_y_var);
                    y_max = std::max(y_max, y_var);
                    y_min = std::min(y_min, y_var);
                }
            }
            fout << ti << " " << y_fit << " " << y_max << " " << y_min << "\n";
        }
    }

    // 3. Log-space data and fit (for a second plot)
    {
        std::ofstream fout("logdata.txt");
        fout << "# t [days]   ln(y)   d(ln(y))\n";
        for (int i = 0; i < n; i++)
            fout << t[i] << " " << lny[i] << " " << dlny[i] << "\n";
    }
    {
        std::ofstream fout("logfit.txt");
        fout << "# t   ln(y)_fit   ln(y)_upper   ln(y)_lower\n";
        for (int i = 0; i <= npts; i++) {
            double ti       = tmin + (tmax - tmin) * i / npts;
            double lny_fit  = fit(ti);

            double ln_max = std::numeric_limits<double>::lowest();
            double ln_min = std::numeric_limits<double>::max();
            for (int s1 = -1; s1 <= 1; s1 += 2) {
                for (int s2 = -1; s2 <= 1; s2 += 2) {
                    double ln_y_var = (fit.c[0] + s1 * fit.dc[0])
                                    + (fit.c[1] + s2 * fit.dc[1]) * ti;
                    ln_max = std::max(ln_max, ln_y_var);
                    ln_min = std::min(ln_min, ln_y_var);
                }
            }
            fout << ti << " " << lny_fit << " " << ln_max << " " << ln_min << "\n";
        }
    }

    std::printf("Data files written: data.txt, fit.txt, logdata.txt, logfit.txt\n");
    return 0;
}
