// "householder.h" header file.
// QR-decomposition via Householder reflections.
// More numerically stable than Gram-Schmidt.
#pragma once

#include "matrix.h"

namespace pp {

class householder_qr {
public:
    pp::matrix Q;
    pp::matrix R;
    int reflection_count; // number of Householder reflections applied

    // Constructor: perform QR-decomposition of A (n×m, n >= m)
    explicit householder_qr(const pp::matrix& A);

    // Solve QRx = b  =>  Rx = Q^T b  (back-substitution)
    pp::vector solve(const pp::vector& b) const;

    // Determinant: det(A) = det(Q) * det(R), det(Q) = (-1)^reflections
    double det() const;

    // Inverse
    pp::matrix inverse() const;
};

} // namespace pp
