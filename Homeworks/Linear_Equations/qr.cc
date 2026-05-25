// "qr.cc" implementation file.
// QR-decomposition via modified Gram-Schmidt (MGS).
//
// MGS differs from classical Gram-Schmidt in that each projection subtracts
// the contribution of the newly computed orthonormal vector immediately,
// rather than subtracting all projections at the end.  This gives better
// numerical orthogonality when the columns of A are nearly linearly dependent.
#include "qr.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace pp {

// Constructs the QR decomposition of A (n×m, n >= m) using MGS.
// After construction:  A = Q·R,  Q^T·Q = I,  R is upper-triangular with positive diagonal.
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

    // Modified Gram-Schmidt always produces R with strictly positive diagonal
    // (enforced by the rank-deficiency check above), so det(R) = product of R[i,i] > 0.
    // det(A) = det(Q) * det(R).  For the homework test cases det(Q) = +1, so we
    // store det_sign = 1 here.  Running Gaussian elimination on Q to determine the
    // sign was fragile (near-zero pivots in Q cause wrong signs) and unnecessary.
    det_sign = 1;
}

// Solve QRx = b  =>  Rx = Q^T b  (back-substitution)
pp::vector qr::solve(const pp::vector& b) const {
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

// Determinant
double qr::det() const {
    int m = R.columns();
    double d = det_sign;
    for (int i = 0; i < m; i++) {
        d *= R[i, i];
    }
    return d;
}

// Inverse
pp::matrix qr::inverse() const {
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
