// "householder.cc" implementation file.
// QR-decomposition via Householder reflections.
//
// Each Householder reflector H_k = I - 2*v*v^T / (v^T*v) zeros out the
// sub-diagonal elements of column k by reflecting x onto alpha*e1, where
// alpha = -sign(x[0]) * ||x||.  Accumulating all reflectors gives the full
// orthogonal factor Q.  Compared to Gram-Schmidt, Householder reflections
// maintain orthogonality to near machine precision even for ill-conditioned A.
#include "householder.h"
#include <cmath>

namespace pp {

// Constructor: perform QR-decomposition of A (n×m, n >= m)
// using Householder reflections.
householder_qr::householder_qr(const pp::matrix& A) {
    int n = A.rows();
    int m = A.columns();

    R = A; // work on a copy
    Q = pp::matrix::identity(n);
    reflection_count = 0;

    for (int k = 0; k < m; k++) {
        // Extract the subcolumn x = R[k:n, k]
        int len = n - k;
        pp::vector x(len);
        for (int i = 0; i < len; i++) {
            x[i] = R[k + i, k];
        }

        // Compute v = x - alpha * e1
        // where alpha = -sign(x[0]) * ||x||
        double x_norm = x.norm();
        if (x_norm < 1e-12) continue; // already zero

        double alpha = (x[0] >= 0.0) ? -x_norm : x_norm;
        pp::vector v(len);
        v[0] = x[0] - alpha;
        for (int i = 1; i < len; i++) {
            v[i] = x[i];
        }

        double v_dot = v.dot(v);
        if (v_dot < 1e-24) continue;

        reflection_count++;

        // Apply Householder reflection to R:
        // R[k:n, k:m] -= 2 * v * (v^T * R[k:n, k:m]) / (v^T v)
        for (int j = k; j < m; j++) {
            double dot = 0.0;
            for (int i = 0; i < len; i++) {
                dot += v[i] * R[k + i, j];
            }
            double coeff = 2.0 * dot / v_dot;
            for (int i = 0; i < len; i++) {
                R[k + i, j] -= coeff * v[i];
            }
        }

        // Accumulate Q: Q[:, k:n] -= 2 * Q[:, k:n] * v * v^T / (v^T v)
        for (int i = 0; i < n; i++) {
            double dot = 0.0;
            for (int j = 0; j < len; j++) {
                dot += Q[i, k + j] * v[j];
            }
            double coeff = 2.0 * dot / v_dot;
            for (int j = 0; j < len; j++) {
                Q[i, k + j] -= coeff * v[j];
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
pp::vector householder_qr::solve(const pp::vector& b) const {
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

// Determinant: det(A) = det(Q) * det(R), det(Q) = (-1)^reflections
double householder_qr::det() const {
    int m = R.columns();
    double d = (reflection_count % 2 == 0) ? 1.0 : -1.0;
    for (int i = 0; i < m; i++) {
        d *= R[i, i];
    }
    return d;
}

// Inverse
pp::matrix householder_qr::inverse() const {
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
