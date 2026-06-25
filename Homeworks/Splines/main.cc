// Nikolaos Antonios Iakynthos Nemet
// This project was created using a combination of Claude Code and GitHub Copilot.

// "main.cc" -- Spline Interpolation.
// Linear, quadratic, and cubic spline interpolation with derivatives
// and integration.  Produces data files for gnuplot visualization.
#include <cstdio>
#include <cmath>
#include <fstream>
#include "vector.h"
#include "spline.h"

int main() {
    std::printf("=== Spline Interpolation ===\n\n");

    int n = 10;  // data points: x = 0, 1, ..., 9

    // -- Task A: Linear Spline -----------------------------------------
    std::printf("-- Task A: Linear Spline ----------------------------------------\n\n");

    // Data: {x_i = 0, 1, ..., 9;  y_i = cos(x_i)}
    pp::vector xcos(n), ycos(n);
    for (int i = 0; i < n; i++) {
        xcos[i] = (NUMBER)i;
        ycos[i] = std::cos((NUMBER)i);
    }
    std::printf("Data: x = 0, 1, ..., 9;  y = cos(x)\n\n");

    // Test interpolation at several points
    NUMBER z1 = 1.5, z2 = 4.7, z3 = 8.3;
    std::printf("  linterp(%.1f) = %+.8f   exact cos(%.1f) = %+.8f\n",
                (double)z1, (double)pp::linterp(xcos, ycos, z1),
                (double)z1, (double)std::cos(z1));
    std::printf("  linterp(%.1f) = %+.8f   exact cos(%.1f) = %+.8f\n",
                (double)z2, (double)pp::linterp(xcos, ycos, z2),
                (double)z2, (double)std::cos(z2));
    std::printf("  linterp(%.1f) = %+.8f   exact cos(%.1f) = %+.8f\n",
                (double)z3, (double)pp::linterp(xcos, ycos, z3),
                (double)z3, (double)std::cos(z3));

    // Test integral: int_0^z cos(x) dx = sin(z)
    NUMBER z4 = 9.0;
    std::printf("\n  linterp_integ(%.1f) = %+.8f   exact sin(%.1f) = %+.8f\n",
                (double)z4, (double)pp::linterp_integ(xcos, ycos, z4),
                (double)z4, (double)std::sin(z4));

    // Write data files
    {
        std::ofstream fout("cos_points.txt");
        fout << "# x   cos(x)\n";
        for (int i = 0; i < n; i++) {
            fout << xcos[i] << " " << ycos[i] << "\n";
        }
    }
    {
        std::ofstream fout("linear_interp.txt");
        fout << "# z   linterp(z)   cos(z)   linterp_integ(z)   sin(z)\n";
        int npts = 500;
        NUMBER zmin = xcos[0], zmax = xcos[n - 1];
        for (int i = 0; i <= npts; i++) {
            NUMBER z = zmin + (zmax - zmin) * i / npts;
            fout << (double)z << " "
                 << (double)pp::linterp(xcos, ycos, z) << " "
                 << (double)std::cos(z) << " "
                 << (double)pp::linterp_integ(xcos, ycos, z) << " "
                 << (double)std::sin(z) << "\n";
        }
    }
    std::printf("\n  Data files: cos_points.txt, linear_interp.txt\n\n");

    // -- Task B: Quadratic Spline --------------------------------------
    std::printf("-- Task B: Quadratic Spline -------------------------------------\n\n");

    // Debug tests with simple polynomial tables (x = 1, 2, 3, 4, 5)
    std::printf("Debug tests (x = 1, 2, 3, 4, 5):\n\n");
    {
        pp::vector xt = {1, 2, 3, 4, 5};

        // y = 1 (constant) -- expect b = {0,...}, c = {0,...}
        pp::vector y1 = {1, 1, 1, 1, 1};
        pp::qspline qs1(xt, y1);
        std::printf("  y = 1:    b = %s\n", qs1.b.to_string().c_str());
        std::printf("            c = %s\n", qs1.c.to_string().c_str());

        // y = x (linear) -- expect b = {1,...}, c = {0,...}
        pp::vector y2 = {1, 2, 3, 4, 5};
        pp::qspline qs2(xt, y2);
        std::printf("  y = x:    b = %s\n", qs2.b.to_string().c_str());
        std::printf("            c = %s\n", qs2.c.to_string().c_str());

        // y = x^2 (quadratic) -- expect b = {2,4,6,8}, c = {1,1,1,1}
        pp::vector y3 = {1, 4, 9, 16, 25};
        pp::qspline qs3(xt, y3);
        std::printf("  y = x^2:  b = %s\n", qs3.b.to_string().c_str());
        std::printf("            c = %s\n", qs3.c.to_string().c_str());
    }
    std::printf("\n");

    // Data: {x_i = 0, 1, ..., 9;  y_i = sin(x_i)}
    pp::vector xsin(n), ysin(n);
    for (int i = 0; i < n; i++) {
        xsin[i] = (NUMBER)i;
        ysin[i] = std::sin((NUMBER)i);
    }

    pp::qspline qs(xsin, ysin);
    std::printf("Data: x = 0, 1, ..., 9;  y = sin(x)\n\n");
    std::printf("  qspline(%.1f)       = %+.8f   exact sin(%.1f) = %+.8f\n",
                (double)z1, (double)qs.evaluate(z1),
                (double)z1, (double)std::sin(z1));
    std::printf("  qspline'(%.1f)      = %+.8f   exact cos(%.1f) = %+.8f\n",
                (double)z1, (double)qs.derivative(z1),
                (double)z1, (double)std::cos(z1));
    // int_0^z sin(x) dx = 1 - cos(z)
    std::printf("  qspline_integ(%.1f) = %+.8f   exact 1-cos(%.1f) = %+.8f\n",
                (double)z4, (double)qs.integral(z4),
                (double)z4, (double)(1.0 - std::cos(z4)));

    // Write data files
    {
        std::ofstream fout("sin_points.txt");
        fout << "# x   sin(x)\n";
        for (int i = 0; i < n; i++) {
            fout << xsin[i] << " " << ysin[i] << "\n";
        }
    }
    {
        std::ofstream fout("quad_interp.txt");
        fout << "# z  qspline(z)  sin(z)  qspline'(z)  cos(z)"
             << "  qspline_integ(z)  1-cos(z)\n";
        int npts = 500;
        NUMBER zmin = xsin[0], zmax = xsin[n - 1];
        for (int i = 0; i <= npts; i++) {
            NUMBER z = zmin + (zmax - zmin) * i / npts;
            fout << (double)z << " "
                 << (double)qs.evaluate(z) << " "
                 << (double)std::sin(z) << " "
                 << (double)qs.derivative(z) << " "
                 << (double)std::cos(z) << " "
                 << (double)qs.integral(z) << " "
                 << (double)(1.0 - std::cos(z)) << "\n";
        }
    }
    std::printf("\n  Data files: sin_points.txt, quad_interp.txt\n\n");

    // -- Task C: Cubic Spline ------------------------------------------
    std::printf("-- Task C: Cubic Spline -----------------------------------------\n\n");

    pp::cspline cs(xcos, ycos);

    std::printf("Data: x = 0, 1, ..., 9;  y = cos(x)\n\n");
    std::printf("  cspline(%.1f)       = %+.8f   exact cos(%.1f)  = %+.8f\n",
                (double)z1, (double)cs.evaluate(z1),
                (double)z1, (double)std::cos(z1));
    std::printf("  cspline'(%.1f)      = %+.8f   exact -sin(%.1f) = %+.8f\n",
                (double)z1, (double)cs.derivative(z1),
                (double)z1, (double)(-std::sin(z1)));
    // int_0^z cos(x) dx = sin(z)
    std::printf("  cspline_integ(%.1f) = %+.8f   exact sin(%.1f)  = %+.8f\n",
                (double)z4, (double)cs.integral(z4),
                (double)z4, (double)std::sin(z4));

    std::printf("\n  Comparison: gnuplot's 'smooth csplines' uses the same natural\n");
    std::printf("  cubic spline -- see plot_cubic.svg for overlay.\n");

    // Write data files for plotting
    {
        std::ofstream fout("cubic_interp.txt");
        fout << "# z  cspline(z)  cos(z)  cspline'(z)  -sin(z)"
             << "  cspline_integ(z)  sin(z)\n";
        int npts = 500;
        NUMBER zmin = xcos[0], zmax = xcos[n - 1];
        for (int i = 0; i <= npts; i++) {
            NUMBER z = zmin + (zmax - zmin) * i / npts;
            fout << (double)z << " "
                 << (double)cs.evaluate(z) << " "
                 << (double)std::cos(z) << " "
                 << (double)cs.derivative(z) << " "
                 << (double)(-std::sin(z)) << " "
                 << (double)cs.integral(z) << " "
                 << (double)std::sin(z) << "\n";
        }
    }
    std::printf("\n  Data files: cubic_interp.txt\n\n");

    // -- Task C (functional style): make_qspline ----------------------
    std::printf("-- Task C: Functional-style quadratic spline (make_qspline) -----\n\n");

    // make_qspline returns a std::function<double(double)> closure that
    // captures all spline data by value.  This is the functional programming
    // alternative to the qspline OOP class: no visible internal state, just
    // a plain callable.
    auto S = pp::make_qspline(xsin, ysin);

    std::printf("  Using make_qspline (functional style, captures x,y,b,c by value):\n\n");
    std::printf("  S(%.1f)  = %+.8f   exact sin(%.1f) = %+.8f\n",
                (double)z1, S(z1), (double)z1, std::sin(z1));
    std::printf("  S(%.1f)  = %+.8f   exact sin(%.1f) = %+.8f\n",
                (double)z2, S(z2), (double)z2, std::sin(z2));

    // Verify the closure gives the same result as the OOP version.
    double oop_val  = qs.evaluate(z1);
    double func_val = S(z1);
    std::printf("\n  Agreement check at z=%.1f: |OOP - closure| = %.2e %s\n\n",
                (double)z1, std::abs(oop_val - func_val),
                std::abs(oop_val - func_val) < 1e-14 ? "OK" : "FAIL");

    std::printf("All data files written successfully.\n");
    return 0;
}
