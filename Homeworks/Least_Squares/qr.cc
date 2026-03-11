// "qr.cc" implementation file.
#include "qr.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace pp {

// -- Constructor: A = QR via modified Gram-Schmidt ----------------
//
// For each column i of A:
//   R_{ii} = ||q_i||           (normalise)
//   q_i   = q_i / R_{ii}
//   R_{ij} = q_i * q_j         (project remaining columns)
//   q_j   = q_j - R_{ij} q_i   (subtract projection)
//
// Modified GS operates on the *partially orthogonalised* columns
// rather than the originals, which improves numerical stability.
//
qr::qr(const pp::matrix& A) {
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

    // Compute sign of det(Q) for square matrices.
    // Gram-Schmidt always produces R_{ii} > 0, so
    //   det(A) = det(Q) * prod(R_{ii})
    // and sign(det A) = sign(det Q).  We find sign(det Q) via
    // an in-place LU with partial pivoting on a copy of Q.
    det_sign = 1;
    if (n == m) {
        pp::matrix W = Q;
        for (int k = 0; k < n; k++) {
            int pivot = k;
            double max_val = std::abs(W[k, k]);
            for (int ii = k + 1; ii < n; ii++) {
                if (std::abs(W[ii, k]) > max_val) {
                    max_val = std::abs(W[ii, k]);
                    pivot = ii;
                }
            }
            if (pivot != k) {
                for (int jj = k; jj < n; jj++)
                    std::swap(W[k, jj], W[pivot, jj]);
                det_sign = -det_sign;
            }
            for (int ii = k + 1; ii < n; ii++) {
                double factor = W[ii, k] / W[k, k];
                for (int jj = k + 1; jj < n; jj++)
                    W[ii, jj] -= factor * W[k, jj];
            }
            if (W[k, k] < 0) det_sign = -det_sign;
        }
    }
}

// -- Solve Ax = b in the least-squares sense --------------------
//
// From A = QR:  Ax = b  =>  QRx = b  =>  Rx = Q^T b.
// Step 1: c = Q^T b  (projection onto the orthonormal basis)
// Step 2: back-substitution on the upper-triangular system Rx = c
//
pp::vector qr::solve(const pp::vector& b) const {
    int m = R.columns();
    // c_i = q_i * b  (Q is stored column-wise)
    pp::vector c(m);
    for (int i = 0; i < m; i++) {
        c[i] = Q[i].dot(b);
    }
    // Back-substitution: x_i = (c_i - sum_{j>i}( R_{ij} x_j )) / R_{ii}
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

// -- Determinant: det(A) = det(Q) * prod_i( R_{ii} ) ----------
double qr::det() const {
    int m = R.columns();
    double d = det_sign;
    for (int i = 0; i < m; i++) {
        d *= R[i, i];
    }
    return d;
}

// -- Inverse: A^{-1} column-by-column via solve(e_j) -----------
pp::matrix qr::inverse() const {
    int m = R.columns();
    pp::matrix B(m, m);
    for (int j = 0; j < m; j++) {
        pp::vector ej(m);
        ej[j] = 1.0;
        B[j] = solve(ej);   // j-th column of A^{-1}
    }
    return B;
}

} // namespace pp
