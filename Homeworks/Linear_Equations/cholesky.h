// "cholesky.h" header file.
// Cholesky decomposition: A = LL^T for symmetric positive definite matrices.
// L is lower triangular.
#pragma once

#include <stdexcept>
#include "matrix.h"

namespace pp {

class cholesky {
public:
    pp::matrix L;

    // Constructor: perform Cholesky decomposition of symmetric positive
    // definite matrix A. Produces lower triangular L such that A = LL^T.
    explicit cholesky(const pp::matrix& A) {
        int n = A.rows();
        if (n != A.columns()) {
            throw std::runtime_error("Cholesky decomposition requires a square matrix");
        }

        L = pp::matrix(n, n);

        for (int j = 0; j < n; j++) {
            // Diagonal element
            double sum = A[j, j];
            for (int k = 0; k < j; k++) {
                sum -= L[j, k] * L[j, k];
            }
            if (sum <= 0.0) {
                throw std::runtime_error("Cholesky decomposition failed: matrix is not positive definite");
            }
            L[j, j] = std::sqrt(sum);

            // Off-diagonal elements in column j
            for (int i = j + 1; i < n; i++) {
                double s = A[i, j];
                for (int k = 0; k < j; k++) {
                    s -= L[i, k] * L[j, k];
                }
                L[i, j] = s / L[j, j];
            }
        }
    }

    // Forward substitution: solve Ly = b
    pp::vector forward_substitute(const pp::vector& b) const {
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

    // Back substitution: solve L^T x = y
    pp::vector back_substitute(const pp::vector& y) const {
        int n = L.rows();
        pp::vector x(n);
        for (int i = n - 1; i >= 0; i--) {
            double sum = y[i];
            for (int j = i + 1; j < n; j++) {
                sum -= L[j, i] * x[j]; // L^T[i,j] = L[j,i]
            }
            x[i] = sum / L[i, i];
        }
        return x;
    }

    // Solve Ax = b  =>  LL^T x = b
    // Forward: Ly = b, Back: L^T x = y
    pp::vector solve(const pp::vector& b) const {
        pp::vector y = forward_substitute(b);
        return back_substitute(y);
    }

    // Determinant: det(A) = det(L)^2 = (product of L diagonal)^2
    double det() const {
        int n = L.rows();
        double d = 1.0;
        for (int i = 0; i < n; i++) {
            d *= L[i, i];
        }
        return d * d;
    }

    // Inverse: solve Ax_j = e_j for each unit vector
    pp::matrix inverse() const {
        int n = L.rows();
        pp::matrix B(n, n);
        for (int j = 0; j < n; j++) {
            pp::vector ej(n);
            ej[j] = 1.0;
            B[j] = solve(ej);
        }
        return B;
    }
};

} // namespace pp
