// "evd.cc" implementation file.
// Jacobi eigenvalue decomposition for real symmetric matrices.
#include "evd.h"
#include <cmath>

namespace pp {

// Multiply A from the right by Jacobi rotation J(p,q,θ): A ← A·J
// The Jacobi matrix J differs from the identity only in the (p,p), (p,q),
// (q,p), (q,q) entries, so only columns p and q of A are modified.
void EVD::timesJ(pp::matrix& A, int p, int q, double theta) {
    double c = std::cos(theta), s = std::sin(theta);
    for (int i = 0; i < A.rows(); i++) {
        double aip = A[i, p], aiq = A[i, q];
        A[i, p] = c * aip - s * aiq;
        A[i, q] = s * aip + c * aiq;
    }
}

// Multiply A from the left by J(p,q,θ)^T: A ← J^T·A
// Only rows p and q of A are modified.
void EVD::Jtimes(pp::matrix& A, int p, int q, double theta) {
    double c = std::cos(theta), s = std::sin(theta);
    for (int j = 0; j < A.columns(); j++) {
        double apj = A[p, j], aqj = A[q, j];
        A[p, j] =  c * apj + s * aqj;
        A[q, j] = -s * apj + c * aqj;
    }
}

// ── Standard Jacobi eigenvalue algorithm with cyclic sweeps ─────
// Convergence criterion: the sweep is repeated until none of the
// diagonal elements change.  The exact floating-point comparison
// (!=) is intentional: atan2 returns θ=0 when apq is negligible
// relative to (aqq−app), making the rotation an exact identity and
// the diagonal elements unchanged bit-for-bit.  Once every (p,q)
// pair produces no change, the matrix is diagonal to machine
// precision and the loop terminates.
EVD::EVD(pp::matrix A) : V(A.rows(), A.rows()), w(A.rows()) {
    int n = A.rows();
    V.setid();

    bool changed;
    do {
        changed = false;
        for (int p = 0; p < n - 1; p++)
            for (int q = p + 1; q < n; q++) {
                double apq = A[p, q], app = A[p, p], aqq = A[q, q];
                double theta = 0.5 * std::atan2(2 * apq, aqq - app);
                double c = std::cos(theta), s = std::sin(theta);
                double new_app = c * c * app - 2 * s * c * apq + s * s * aqq;
                double new_aqq = s * s * app + 2 * s * c * apq + c * c * aqq;
                if (new_app != app || new_aqq != aqq) {
                    changed = true;
                    timesJ(A, p, q,  theta); // A ← A·J
                    Jtimes(A, p, q, -theta); // A ← J^T·A
                    timesJ(V, p, q,  theta); // V ← V·J
                }
            }
    } while (changed);

    // Copy diagonal elements into w
    for (int i = 0; i < n; i++) w[i] = A[i, i];
}

// ── Optimized Jacobi: update only upper triangle of A ───────────
// Since A is symmetric we only need to maintain the upper triangle.
// The rotation formulae for the upper-triangle elements are:
//   For i ≠ p,q:  A'[i,p] = c·A[i,p] − s·A[i,q]
//                 A'[i,q] = s·A[i,p] + c·A[i,q]
// applied only to the stored upper-triangle positions.  The diagonal
// elements are updated analytically.  This saves roughly half the
// work compared to the full timesJ/Jtimes approach.
EVD EVD::optimized(pp::matrix A) {
    return EVD(std::move(A), OptimizedTag{});
}

EVD::EVD(pp::matrix A, OptimizedTag) : V(A.rows(), A.rows()), w(A.rows()) {
    int n = A.rows();
    V.setid();

    bool changed;
    do {
        changed = false;
        for (int p = 0; p < n - 1; p++)
            for (int q = p + 1; q < n; q++) {
                double apq = A[p, q], app = A[p, p], aqq = A[q, q];
                double theta = 0.5 * std::atan2(2 * apq, aqq - app);
                double c = std::cos(theta), s = std::sin(theta);
                double new_app = c * c * app - 2 * s * c * apq + s * s * aqq;
                double new_aqq = s * s * app + 2 * s * c * apq + c * c * aqq;
                if (new_app != app || new_aqq != aqq) {
                    changed = true;
                    // Update only the upper triangle of A (symmetric)
                    // The full transform is A' = J^T·A·J.  For i ≠ p,q:
                    //   A'[i,p] = c·A[i,p] − s·A[i,q]
                    //   A'[i,q] = s·A[i,p] + c·A[i,q]
                    // We map these to upper-triangle storage positions.

                    // For i < p: (i,p) and (i,q) are both upper-triangle
                    for (int i = 0; i < p; i++) {
                        double aip = A[i, p], aiq = A[i, q];
                        A[i, p] = c * aip - s * aiq;
                        A[i, q] = s * aip + c * aiq;
                    }
                    // For p < i < q: A[i,p]=A[p,i] (upper), A[i,q] (upper)
                    for (int i = p + 1; i < q; i++) {
                        double api = A[p, i], aiq = A[i, q]; // api = A[i,p] by symmetry
                        A[p, i] = c * api - s * aiq;
                        A[i, q] = s * api + c * aiq;
                    }
                    // For i > q: A[i,p]=A[p,i] (upper), A[i,q]=A[q,i] (upper)
                    for (int i = q + 1; i < n; i++) {
                        double api = A[p, i], aqi = A[q, i]; // api=A[i,p], aqi=A[i,q]
                        A[p, i] = c * api - s * aqi;
                        A[q, i] = s * api + c * aqi;
                    }
                    // Update diagonal elements analytically
                    A[p, p] = new_app;
                    A[q, q] = new_aqq;
                    A[p, q] = 0.0; // the rotation zeroes out this element
                    // Update V (full matrix, not symmetric)
                    timesJ(V, p, q, theta);
                }
            }
    } while (changed);

    for (int i = 0; i < n; i++) w[i] = A[i, i];
}

} // namespace pp
