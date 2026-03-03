// "qr.h" header file. QR-decomposition via modified Gram-Schmidt.
#pragma once

#include <stdexcept>
#include "matrix.h"

namespace pp {

class qr {
public:
    pp::matrix Q;
    pp::matrix R;

    // Constructor: perform QR-decomposition of A (n×m, n >= m)
    // using modified Gram-Schmidt orthogonalization.
    explicit qr(const pp::matrix& A) {
        int n = A.rows();
        int m = A.columns();
        Q = pp::matrix(n, m);
        R = pp::matrix(m, m);

        // Copy columns of A into Q
        for (int j = 0; j < m; j++) {
            Q[j] = A[j];
        }

        // Modified Gram-Schmidt
        for (int i = 0; i < m; i++) {
            R[i, i] = Q[i].norm();
            if (R[i, i] < 1e-12) {
                throw std::runtime_error("QR decomposition failed: matrix is rank-deficient");
            }
            Q[i] = Q[i] / R[i, i];
            for (int j = i + 1; j < m; j++) {
                R[i, j] = Q[i].dot(Q[j]);
                Q[j] = Q[j] - Q[i] * R[i, j];
            }
        }
    }

    // Solve QRx = b  =>  Rx = Q^T b  (back-substitution)
    pp::vector solve(const pp::vector& b) const {
        int m = R.columns();
        // Compute c = Q^T * b
        pp::vector c(m);
        for (int i = 0; i < m; i++) {
            c[i] = Q[i].dot(b);
        }
        // Back-substitution on Rx = c
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

    // Determinant of the original matrix (product of R's diagonal)
    double det() const {
        int m = R.columns();
        double d = 1.0;
        for (int i = 0; i < m; i++) {
            d *= R[i, i];
        }
        return d;
    }

    // Inverse of the original matrix A:
    // A = QR  =>  A^{-1} = R^{-1} Q^T
    // Computed by solving A x_j = e_j for each unit vector e_j.
    pp::matrix inverse() const {
        int m = R.columns();
        pp::matrix B(m, m);
        for (int j = 0; j < m; j++) {
            pp::vector ej(m);
            ej[j] = 1.0;
            B[j] = solve(ej);
        }
        return B;
    }
};

} // namespace pp
