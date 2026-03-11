// "qr.h" header file. QR-decomposition via modified Gram-Schmidt.
#pragma once

#include "matrix.h"

namespace pp {

// -----------------------------------------------------------------
// pp::qr -- thin QR-decomposition via modified Gram-Schmidt
//
// Given an nxm matrix A (n >= m, full column rank), the constructor
// computes A = Q R, where
//   Q  is nxm with orthonormal columns (Q^T Q = I_m),
//   R  is mxm upper-triangular with positive diagonal.
//
// Modified Gram-Schmidt is used for better numerical stability
// compared to classical Gram-Schmidt: each column is orthogonalised
// against the *already-updated* previous columns rather than
// against the originals.
//
// Provided operations for the decomposed system:
//   solve(b)  -- least-squares solution of Ax ~ b  (via R x = Q^T b)
//   det()     -- determinant of A  (= det_sign * prod(R_ii))
//   inverse() -- A^{-1} for square A  (column-by-column back-sub)
// -----------------------------------------------------------------
class qr {
public:
    pp::matrix Q;
    pp::matrix R;
    int det_sign; // sign of det(Q): +1 or -1

    // Constructor: perform QR-decomposition of A (nxm, n >= m)
    explicit qr(const pp::matrix& A);

    // Solve QRx = b  =>  Rx = Q^T b  (back-substitution)
    pp::vector solve(const pp::vector& b) const;

    // Determinant of the original matrix: det(Q) * product(R diagonal)
    double det() const;

    // Inverse of the original matrix
    pp::matrix inverse() const;
};

} // namespace pp
