// "lsfit.cc" implementation file.
// Ordinary least-squares fit by QR-decomposition.
#include "lsfit.h"
#include "qr.h"
#include <cmath>
#include <stdexcept>

namespace pp {

// -- Constructor: build the design matrix and solve via QR --------
//
// Given data {x_i, y_i, dy_i} and basis functions {f_k}, form:
//   A_{ik} = f_k(x_i) / dy_i      (nxm design matrix)
//   b_i    = y_i / dy_i            (weighted data vector)
//
// The least-squares solution minimises ||Ac - b||^2 and is obtained
// by QR-decomposing the tall matrix A = QR, then solving Rc = Q^T b.
//
// The covariance matrix of the coefficients is  Cov = (A^T A)^{-1}
// = (R^T R)^{-1} = R^{-1} (R^{-1})^T .
//
lsfit::lsfit(
    const BasisFunctions& fs,
    const pp::vector& x,
    const pp::vector& y,
    const pp::vector& dy)
    : fs(fs)
{
    int n = x.size();                       // number of data points
    int m = static_cast<int>(this->fs.size());  // number of basis functions

    if (n < m)
        throw std::runtime_error("lsfit: more basis functions than data points");

    // Build the design matrix A (nxm) and weighted data vector b (n)
    pp::matrix A(n, m);
    pp::vector b(n);
    for (int i = 0; i < n; i++) {
        b[i] = y[i] / dy[i];
        for (int k = 0; k < m; k++) {
            A[i, k] = this->fs[k](x[i]) / dy[i];
        }
    }

    // QR-decompose the design matrix and solve for coefficients
    pp::qr qr(A);
    c = qr.solve(b);

    // Covariance matrix: Cov = (A^T A)^{-1} = R^{-1} (R^{-1})^T
    //
    // We need R^{-1}, not the pseudo-inverse A^+ = R^{-1} Q^T.
    // Since qr.solve(b) computes R^{-1} Q^T b (the full least-squares
    // solution), it cannot be used directly here.  Instead we perform
    // back-substitution on R alone: solve  R x = e_j  for each unit
    // vector to obtain the columns of R^{-1}.
    pp::matrix Rinv(m, m);
    for (int j = 0; j < m; j++) {
        pp::vector ej(m);
        ej[j] = 1.0;
        pp::vector col(m);
        for (int i = m - 1; i >= 0; i--) {
            double sum = ej[i];
            for (int k = i + 1; k < m; k++)
                sum -= qr.R[i, k] * col[k];
            col[i] = sum / qr.R[i, i];
        }
        Rinv[j] = col;
    }
    cov = Rinv * Rinv.T();

    // Uncertainties: dc_k = sqrt(Cov_{kk})
    dc = pp::vector(m);
    for (int k = 0; k < m; k++)
        dc[k] = std::sqrt(cov[k, k]);
}

// -- Evaluate the fit at a point ----------------------------------
double lsfit::operator()(double x) const {
    double sum = 0.0;
    for (int k = 0; k < c.size(); k++)
        sum += c[k] * fs[k](x);
    return sum;
}

} // namespace pp
