// "lu.cc" implementation file.
#include "lu.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace pp {

// Constructor: perform LU-decomposition of square matrix A
// using Doolittle algorithm with partial pivoting.
lu::lu(const pp::matrix& A) {
    int n = A.rows();
    if (n != A.columns()) {
        throw std::runtime_error("LU decomposition requires a square matrix");
    }

    // Work on a copy
    pp::matrix W = A;
    L = pp::matrix(n, n);
    perm.resize(n);
    sign = 1;

    // Initialize permutation to identity
    for (int i = 0; i < n; i++) perm[i] = i;

    // Initialize L diagonal to 1
    for (int i = 0; i < n; i++) L[i, i] = 1.0;

    // Doolittle with partial pivoting
    for (int k = 0; k < n; k++) {
        // Find pivot: row with largest |W[i,k]| for i = k..n-1
        int pivot = k;
        double max_val = std::abs(W[k, k]);
        for (int i = k + 1; i < n; i++) {
            double val = std::abs(W[i, k]);
            if (val > max_val) {
                max_val = val;
                pivot = i;
            }
        }

        if (max_val < 1e-12) {
            throw std::runtime_error("LU decomposition failed: matrix is singular");
        }

        // Swap rows k and pivot in W, L, and permutation
        if (pivot != k) {
            std::swap(perm[k], perm[pivot]);
            sign = -sign;

            // Swap rows in W (swap elements in each column)
            for (int j = 0; j < n; j++) {
                std::swap(W[k, j], W[pivot, j]);
            }

            // Swap already-computed L entries (columns 0..k-1)
            for (int j = 0; j < k; j++) {
                std::swap(L[k, j], L[pivot, j]);
            }
        }

        // Eliminate below diagonal
        for (int i = k + 1; i < n; i++) {
            L[i, k] = W[i, k] / W[k, k];
            for (int j = k; j < n; j++) {
                W[i, j] -= L[i, k] * W[k, j];
            }
        }
    }

    U = W; // W now contains U
}

// Forward substitution: solve Ly = b
pp::vector lu::forward_substitute(const pp::vector& b) const {
    int n = L.rows();
    pp::vector y(n);
    for (int i = 0; i < n; i++) {
        double sum = b[i];
        for (int j = 0; j < i; j++) {
            sum -= L[i, j] * y[j];
        }
        y[i] = sum / L[i, i];
    }
    return y;
}

// Back substitution: solve Ux = y
pp::vector lu::back_substitute(const pp::vector& y) const {
    int n = U.rows();
    pp::vector x(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = y[i];
        for (int j = i + 1; j < n; j++) {
            sum -= U[i, j] * x[j];
        }
        x[i] = sum / U[i, i];
    }
    return x;
}

// Solve Ax = b  =>  PAx = Pb  =>  LUx = Pb
pp::vector lu::solve(const pp::vector& b) const {
    int n = L.rows();
    pp::vector pb(n);
    for (int i = 0; i < n; i++) {
        pb[i] = b[perm[i]];
    }
    pp::vector y = forward_substitute(pb);
    return back_substitute(y);
}

// Determinant: det(A) = sign * product of U diagonal
double lu::det() const {
    int n = U.rows();
    double d = sign;
    for (int i = 0; i < n; i++) {
        d *= U[i, i];
    }
    return d;
}

// Inverse: solve Ax_j = e_j for each unit vector
pp::matrix lu::inverse() const {
    int n = U.rows();
    pp::matrix B(n, n);
    for (int j = 0; j < n; j++) {
        pp::vector ej(n);
        ej[j] = 1.0;
        B[j] = solve(ej);
    }
    return B;
}

} // namespace pp
