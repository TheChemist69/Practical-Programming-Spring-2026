// "qr.cc" implementation file.
// Modified Gram-Schmidt QR for dense square systems.

#include "qr.h"

#include <cmath>

namespace pp {

QRDecomposition qr_decompose(const matrix& A) {
    const int m = A.rows();
    const int n = A.cols();

    QRDecomposition out;
    out.Q = matrix(m, n);
    out.R = matrix(n, n);
    out.success = true;

    matrix work = A;

    for (int k = 0; k < n; k++) {
        vector ak = work.column(k);
        const double rkk = ak.norm();

        out.R(k, k) = rkk;
        if (rkk <= 1e-14) {
            out.success = false;
            continue;
        }

        const vector qk = ak / rkk;
        out.Q.set_column(k, qk);

        for (int j = k + 1; j < n; j++) {
            vector aj = work.column(j);
            const double rkj = dot(qk, aj);
            out.R(k, j) = rkj;
            aj -= qk * rkj;
            work.set_column(j, aj);
        }
    }

    return out;
}

vector qr_solve(const QRDecomposition& qr, const vector& b, bool* success) {
    const int n = qr.R.rows();
    vector x(n);

    if (!qr.success || qr.Q.rows() != b.size()) {
        if (success) *success = false;
        return x;
    }

    // Compute y = Q^T b.
    vector y(n);
    for (int i = 0; i < n; i++) {
        const vector qi = qr.Q.column(i);
        y[i] = dot(qi, b);
    }

    // Back-substitution in R x = y.
    for (int i = n - 1; i >= 0; i--) {
        double rhs = y[i];
        for (int j = i + 1; j < n; j++) {
            rhs -= qr.R(i, j) * x[j];
        }
        const double rii = qr.R(i, i);
        if (std::abs(rii) <= 1e-14) {
            if (success) *success = false;
            return vector(n);
        }
        x[i] = rhs / rii;
    }

    if (success) *success = true;
    return x;
}

vector solve_linear_system_qr(const matrix& A, const vector& b, bool* success) {
    const QRDecomposition qr = qr_decompose(A);
    return qr_solve(qr, b, success);
}

} // namespace pp
