// "givens.h" header file.
// QR-decomposition via Givens rotations.
// Zeros out sub-diagonal elements one at a time using plane rotations.
#pragma once

#include "matrix.h"

namespace pp {

class givens_qr {
public:
    pp::matrix Q;
    pp::matrix R;

    // Constructor: perform QR-decomposition of A (n×m, n >= m)
    explicit givens_qr(const pp::matrix& A);

    // Solve QRx = b  =>  Rx = Q^T b  (back-substitution)
    pp::vector solve(const pp::vector& b) const;

    // Determinant (product of R diagonal)
    double det() const;

    // Inverse
    pp::matrix inverse() const;
};

} // namespace pp
