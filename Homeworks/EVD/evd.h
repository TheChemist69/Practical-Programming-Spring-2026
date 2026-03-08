// "evd.h" header file.
// Jacobi eigenvalue decomposition for real symmetric matrices.
#pragma once

#include "matrix.h"
#include "vector.h"

namespace pp {

struct EVD {
    pp::matrix V;   // eigenvectors as columns
    pp::vector w;   // eigenvalues

    // Multiply A from the right by Jacobi rotation J(p,q,theta): A <- A*J
    static void timesJ(pp::matrix& A, int p, int q, double theta);

    // Multiply A from the left by Jacobi rotation J(p,q,theta): A <- J^T*A
    static void Jtimes(pp::matrix& A, int p, int q, double theta);

    // Constructor: performs Jacobi diagonalization of symmetric matrix A
    EVD(pp::matrix A);
};

} // namespace pp
