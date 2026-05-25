// "cholesky.h" header file.
// Cholesky decomposition: A = LL^T for symmetric positive definite matrices.
// L is lower triangular.
#pragma once

#include "matrix.h"

namespace pp {

class cholesky {
public:
    pp::matrix L;

    // Constructor: perform Cholesky decomposition of symmetric positive
    // definite matrix A.
    explicit cholesky(const pp::matrix& A);

    // Forward substitution: solve Ly = b
    pp::vector forward_substitute(const pp::vector& b) const;

    // Back substitution: solve L^T x = y
    pp::vector back_substitute(const pp::vector& y) const;

    // Solve Ax = b  =>  LL^T x = b
    pp::vector solve(const pp::vector& b) const;

    // Determinant: det(A) = (product of L diagonal)^2
    double det() const;

    // Inverse
    pp::matrix inverse() const;
};

} // namespace pp
