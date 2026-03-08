// "main.cc" — Jacobi eigenvalue decomposition: tests and hydrogen atom.
#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <fstream>
#include <string>
#include "vector.h"
#include "matrix.h"
#include "evd.h"

// ──────────────────────────────────────────────────────────────────
// Task A: Test the Jacobi diagonalization on a random symmetric matrix
// ──────────────────────────────────────────────────────────────────
void task_A() {
    printf("=== Task A: Jacobi diagonalization test ===\n\n");

    int n = 5;
    // Generate a random symmetric matrix
    pp::matrix R = pp::matrix::random(n, n);
    pp::matrix A = R + R.T(); // symmetric
    printf("Random symmetric matrix A (%dx%d):\n", n, n);
    A.print("");

    // Diagonalize
    pp::EVD evd(A);

    printf("Eigenvalues w:\n");
    evd.w.print("  ");

    // Build diagonal matrix D from eigenvalues
    pp::matrix D(n, n);
    for (int i = 0; i < n; i++) D[i, i] = evd.w[i];

    // Check V^T * A * V == D
    pp::matrix VtAV = evd.V.T() * A * evd.V;
    printf("V^T * A * V (should be D):\n");
    VtAV.print("");

    printf("D:\n");
    D.print("");

    pp::matrix diff1 = VtAV - D;
    printf("||V^T*A*V - D|| = %g\n\n", diff1.norm());

    // Check V * D * V^T == A
    pp::matrix VDVT = evd.V * D * evd.V.T();
    printf("V*D*V^T (should be A):\n");
    VDVT.print("");

    pp::matrix diff2 = VDVT - A;
    printf("||V*D*V^T - A|| = %g\n\n", diff2.norm());

    // Check V^T * V == I
    pp::matrix VtV = evd.V.T() * evd.V;
    pp::matrix I = pp::matrix::identity(n);
    pp::matrix diff3 = VtV - I;
    printf("V^T*V (should be I):\n");
    VtV.print("");
    printf("||V^T*V - I|| = %g\n\n", diff3.norm());

    // Check V * V^T == I
    pp::matrix VVt = evd.V * evd.V.T();
    pp::matrix diff4 = VVt - I;
    printf("V*V^T (should be I):\n");
    VVt.print("");
    printf("||V*V^T - I|| = %g\n\n", diff4.norm());
}

// ──────────────────────────────────────────────────────────────────
// Task B: Hydrogen atom s-wave radial Schrödinger equation
// ──────────────────────────────────────────────────────────────────
void task_B(double rmax, double dr) {
    printf("=== Task B: Hydrogen atom (rmax=%.2f, dr=%.4f) ===\n\n", rmax, dr);

    int npoints = (int)(rmax / dr) - 1;
    if (npoints < 2) {
        printf("Too few grid points (%d). Increase rmax or decrease dr.\n", npoints);
        return;
    }

    // Build the grid
    pp::vector r(npoints);
    for (int i = 0; i < npoints; i++) r[i] = dr * (i + 1);

    // Build the Hamiltonian matrix H = K + W
    pp::matrix H(npoints, npoints);
    for (int i = 0; i < npoints - 1; i++) {
        H[i, i]     = -2 * (-0.5 / dr / dr);
        H[i, i + 1] =  1 * (-0.5 / dr / dr);
        H[i + 1, i] =  1 * (-0.5 / dr / dr);
    }
    H[npoints - 1, npoints - 1] = -2 * (-0.5 / dr / dr);
    for (int i = 0; i < npoints; i++) H[i, i] += -1.0 / r[i];

    // Diagonalize
    pp::EVD evd(H);

    // Print lowest eigenvalues and compare with exact: E_k = -1/(2k^2)
    int nshow = std::min(npoints, 10);
    // Sort eigenvalues to find lowest
    // EVD eigenvalues may not be sorted; find the indices of the smallest ones
    std::vector<int> idx(npoints);
    for (int i = 0; i < npoints; i++) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
        return evd.w[a] < evd.w[b];
    });

    printf("Lowest eigenvalues (k=1..%d):\n", nshow);
    printf("%5s %15s %15s %15s\n", "k", "numerical", "exact", "difference");
    for (int k = 0; k < nshow; k++) {
        double exact = -0.5 / ((k + 1.0) * (k + 1.0));
        double num = evd.w[idx[k]];
        printf("%5d %15.8f %15.8f %15.2e\n", k + 1, num, exact, num - exact);
    }
    printf("\n");

    // Write wave functions to file for plotting
    {
        std::ofstream fout("wavefunctions.data");
        fout << "# r  f1(r)  f2(r)  f3(r)  exact_f1(r)  exact_f2(r)  exact_f3(r)\n";
        int nwf = std::min(3, npoints);
        double Const = 1.0 / std::sqrt(dr);
        for (int i = 0; i < npoints; i++) {
            fout << r[i];
            for (int k = 0; k < nwf; k++) {
                // Ensure sign convention: first non-negligible component positive
                double sign = 1.0;
                for (int ii = 0; ii < npoints; ii++) {
                    if (std::abs(evd.V[ii, idx[k]]) > 1e-6) {
                        sign = (evd.V[ii, idx[k]] > 0) ? 1.0 : -1.0;
                        break;
                    }
                }
                fout << " " << sign * Const * evd.V[i, idx[k]];
            }
            // Exact reduced radial wave functions f_k(r) = r * R_{k+1,0}(r)
            // For hydrogen: f_n(r) = r * R_{n0}(r)
            // n=1: f(r) = 2*r*exp(-r)
            // n=2: f(r) = r/(2*sqrt(2)) * (2 - r) * exp(-r/2)
            // n=3: f(r) = 2*r/(81*sqrt(3)) * (27 - 18*r + 2*r^2) * exp(-r/3)
            double ri = r[i];
            double f1_exact = 2.0 * ri * std::exp(-ri);
            double f2_exact = ri / (2.0 * std::sqrt(2.0)) * (2.0 - ri) * std::exp(-ri / 2.0);
            double f3_exact = 2.0 * ri / (81.0 * std::sqrt(3.0)) * (27.0 - 18.0 * ri + 2.0 * ri * ri) * std::exp(-ri / 3.0);
            fout << " " << f1_exact << " " << f2_exact << " " << f3_exact;
            fout << "\n";
        }
    }
    printf("Wave function data written to wavefunctions.data\n\n");
}

// ──────────────────────────────────────────────────────────────────
// Convergence study: vary dr at fixed rmax
// ──────────────────────────────────────────────────────────────────
void convergence_dr(double rmax) {
    std::ofstream fout("convergence_dr.data");
    fout << "# dr  E0_numerical  E0_exact  error\n";
    double exact = -0.5;
    for (double dr = 0.05; dr <= 1.01; dr += 0.05) {
        int npoints = (int)(rmax / dr) - 1;
        if (npoints < 2) continue;
        pp::vector r(npoints);
        for (int i = 0; i < npoints; i++) r[i] = dr * (i + 1);
        pp::matrix H(npoints, npoints);
        for (int i = 0; i < npoints - 1; i++) {
            H[i, i]     = -2 * (-0.5 / dr / dr);
            H[i, i + 1] =  1 * (-0.5 / dr / dr);
            H[i + 1, i] =  1 * (-0.5 / dr / dr);
        }
        H[npoints - 1, npoints - 1] = -2 * (-0.5 / dr / dr);
        for (int i = 0; i < npoints; i++) H[i, i] += -1.0 / r[i];
        pp::EVD evd(H);
        // Find lowest eigenvalue
        double e0 = evd.w[0];
        for (int i = 1; i < npoints; i++) if (evd.w[i] < e0) e0 = evd.w[i];
        fout << dr << " " << e0 << " " << exact << " " << std::abs(e0 - exact) << "\n";
    }
    printf("Convergence (dr) data written to convergence_dr.data\n");
}

// ──────────────────────────────────────────────────────────────────
// Convergence study: vary rmax at fixed dr
// ──────────────────────────────────────────────────────────────────
void convergence_rmax(double dr) {
    std::ofstream fout("convergence_rmax.data");
    fout << "# rmax  E0_numerical  E0_exact  error\n";
    double exact = -0.5;
    for (double rmax = 1.0; rmax <= 20.01; rmax += 1.0) {
        int npoints = (int)(rmax / dr) - 1;
        if (npoints < 2) continue;
        pp::vector r(npoints);
        for (int i = 0; i < npoints; i++) r[i] = dr * (i + 1);
        pp::matrix H(npoints, npoints);
        for (int i = 0; i < npoints - 1; i++) {
            H[i, i]     = -2 * (-0.5 / dr / dr);
            H[i, i + 1] =  1 * (-0.5 / dr / dr);
            H[i + 1, i] =  1 * (-0.5 / dr / dr);
        }
        H[npoints - 1, npoints - 1] = -2 * (-0.5 / dr / dr);
        for (int i = 0; i < npoints; i++) H[i, i] += -1.0 / r[i];
        pp::EVD evd(H);
        double e0 = evd.w[0];
        for (int i = 1; i < npoints; i++) if (evd.w[i] < e0) e0 = evd.w[i];
        fout << rmax << " " << e0 << " " << exact << " " << std::abs(e0 - exact) << "\n";
    }
    printf("Convergence (rmax) data written to convergence_rmax.data\n");
}

// ──────────────────────────────────────────────────────────────────
// Timing / scaling study (Task C)
// ──────────────────────────────────────────────────────────────────
void task_C_timing(int N) {
    // Diagonalize a random symmetric NxN matrix; used for timing
    pp::matrix R = pp::matrix::random(N, N);
    pp::matrix A = R + R.T();
    pp::EVD evd(A);
    // Print the first eigenvalue so the compiler doesn't optimize away
    printf("%d %g\n", N, evd.w[0]);
}

// ──────────────────────────────────────────────────────────────────
int main(int argc, char* argv[]) {
    // Parse command line
    double rmax = 10.0;
    double dr = 0.3;
    std::string mode = "all";

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "-rmax") == 0 && i + 1 < argc) {
            rmax = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "-dr") == 0 && i + 1 < argc) {
            dr = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "-mode") == 0 && i + 1 < argc) {
            mode = argv[++i];
        } else if (std::strcmp(argv[i], "-size") == 0 && i + 1 < argc) {
            int N = std::atoi(argv[++i]);
            task_C_timing(N);
            return 0;
        }
    }

    if (mode == "all" || mode == "taskA") {
        task_A();
    }

    if (mode == "all" || mode == "taskB") {
        task_B(rmax, dr);
    }

    if (mode == "all" || mode == "convergence") {
        convergence_dr(rmax);
        convergence_rmax(dr);
    }

    if (mode == "convergence_dr") {
        convergence_dr(rmax);
    }

    if (mode == "convergence_rmax") {
        convergence_rmax(dr);
    }

    return 0;
}
