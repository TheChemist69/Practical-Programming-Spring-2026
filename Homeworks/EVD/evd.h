// "evd.h" header file.
// Jacobi eigenvalue decomposition for real symmetric matrices.
#pragma once

#include "matrix.h"
#include "vector.h"

namespace pp {

struct EVD {
    pp::matrix V;   // eigenvectors as columns
    pp::vector w;   // eigenvalues

    // Multiply A from the right by Jacobi rotation J(p,q,θ): A ← A·J
    // Only columns p and q are affected → O(n) operations.
    static void timesJ(pp::matrix& A, int p, int q, double theta);

    // Multiply A from the left by J(p,q,θ)^T: A ← J^T·A
    // Only rows p and q are affected → O(n) operations.
    static void Jtimes(pp::matrix& A, int p, int q, double theta);

    // Constructor: performs Jacobi diagonalization of symmetric matrix A
    // using cyclic sweeps over all upper-triangle pairs (p,q).
    EVD(pp::matrix A);

    // Optimized constructor: updates only the upper triangle during
    // rotations, giving roughly 2× speedup for large matrices.
    // Use EVD::optimized(A) factory to invoke.
    static EVD optimized(pp::matrix A);

private:
    // Private tag constructor used by the optimized factory
    struct OptimizedTag {};
    EVD(pp::matrix A, OptimizedTag);
};

} // namespace pp
