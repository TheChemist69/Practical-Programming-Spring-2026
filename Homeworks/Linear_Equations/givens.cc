// "givens.cc" implementation file.
// QR-decomposition via Givens rotations.
//
// Each Givens rotation G(i-1, i, theta) is a plane rotation in the (i-1, i)
// coordinate plane chosen so that R[i, j] becomes zero.  Rotations are applied
// column-by-column, bottom-to-top, zeroing one sub-diagonal element at a time.
// The rotation parameters are c = a/r, s = b/r with r = hypot(a, b), which
// avoids overflow and keeps c^2 + s^2 = 1 numerically.
#include "givens.h"
#include <cmath>

namespace pp {

// Constructor: perform QR-decomposition of A (n×m, n >= m)
// using Givens rotations.
givens_qr::givens_qr(const pp::matrix& A) {
    int n = A.rows();
    int m = A.columns();

    R = A; // work on a copy
    Q = pp::matrix::identity(n);

    for (int j = 0; j < m; j++) {
        // Zero out elements below the diagonal in column j
        for (int i = n - 1; i > j; i--) {
            double a = R[i - 1, j];
            double b = R[i, j];

            if (std::abs(b) < 1e-15) continue; // already zero

            // Compute Givens rotation parameters
            double r = std::hypot(a, b);
            double c = a / r;
            double s = b / r;

            // Apply rotation to R: rows i-1 and i, columns j..m-1
            for (int k = j; k < m; k++) {
                double r1 = R[i - 1, k];
                double r2 = R[i, k];
                R[i - 1, k] =  c * r1 + s * r2;
                R[i, k]     = -s * r1 + c * r2;
            }

            // Accumulate Q: columns i-1 and i
            for (int k = 0; k < n; k++) {
                double q1 = Q[k, i - 1];
                double q2 = Q[k, i];
                Q[k, i - 1] =  c * q1 + s * q2;
                Q[k, i]     = -s * q1 + c * q2;
            }
        }
    }

    // Trim Q to n×m and R to m×m for thin QR
    pp::matrix Q_thin(n, m);
    pp::matrix R_thin(m, m);
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            Q_thin[i, j] = Q[i, j];
        }
        for (int i = 0; i < m; i++) {
            R_thin[i, j] = R[i, j];
        }
    }
    Q = Q_thin;
    R = R_thin;
}

// Solve QRx = b  =>  Rx = Q^T b  (back-substitution)
pp::vector givens_qr::solve(const pp::vector& b) const {
    int m = R.columns();
    pp::vector c(m);
    for (int i = 0; i < m; i++) {
        c[i] = Q[i].dot(b);
    }
    pp::vector x(m);
    for (int i = m - 1; i >= 0; i--) {
        double sum = c[i];
        for (int j = i + 1; j < m; j++) {
            sum -= R[i, j] * x[j];
        }
        x[i] = sum / R[i, i];
    }
    return x;
}

// Determinant (product of R diagonal)
double givens_qr::det() const {
    int m = R.columns();
    double d = 1.0;
    for (int i = 0; i < m; i++) {
        d *= R[i, i];
    }
    return d;
}

// Inverse
pp::matrix givens_qr::inverse() const {
    int m = R.columns();
    pp::matrix B(m, m);
    for (int j = 0; j < m; j++) {
        pp::vector ej(m);
        ej[j] = 1.0;
        B[j] = solve(ej);
    }
    return B;
}

} // namespace pp
