// "qr.h" header file.
// QR-decomposition via modified Gram-Schmidt orthogonalization.
//
// For an n×m matrix A (n >= m), the decomposition A = Q·R produces:
//   Q: n×m matrix with orthonormal columns  (Q^T Q = I_m)
//   R: m×m upper-triangular matrix with positive diagonal
//
// Modified Gram-Schmidt is numerically more stable than classical GS
// because it orthogonalizes each column against the already-processed set
// rather than against all original columns simultaneously.
//
// Applications exposed here:
//   solve(b)   -- least-squares / exact solution of Ax = b
//   det()      -- determinant via product of R diagonal entries
//   inverse()  -- matrix inverse via repeated solve with unit vectors
#pragma once

#include "matrix.h"

namespace pp {

class qr {
public:
    pp::matrix Q;     // orthonormal columns (Q^T Q = I)
    pp::matrix R;     // upper-triangular factor with positive diagonal
    int det_sign;     // +1 always (Gram-Schmidt gives positive R diagonal, so det(A) = product of R[i,i])

    // Constructs the QR decomposition of A (n×m, n >= m).
    // Throws if A is rank-deficient (any R diagonal < 1e-12).
    explicit qr(const pp::matrix& A);

    // Solves QR·x = b by computing c = Q^T·b then back-substituting R·x = c.
    // For square A this gives the exact solution; for tall A it gives the
    // least-squares solution minimising ||Ax - b||.
    pp::vector solve(const pp::vector& b) const;

    // Returns det(A) = det_sign * product(R[i,i]).
    // Since modified Gram-Schmidt produces positive R diagonal, this equals
    // det(R) = det(A) / det(Q).  For most homework matrices det(Q) = +1.
    double det() const;

    // Returns the inverse of A by solving A·x_j = e_j for each unit vector.
    // Only meaningful for square, non-singular A.
    pp::matrix inverse() const;
};

} // namespace pp
