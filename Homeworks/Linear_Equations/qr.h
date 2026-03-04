// "qr.h" header file. QR-decomposition via modified Gram-Schmidt.
#pragma once

#include "matrix.h"

namespace pp {

class qr {
public:
    pp::matrix Q;
    pp::matrix R;
    int det_sign; // sign of det(Q): +1 or -1

    // Constructor: perform QR-decomposition of A (n×m, n >= m)
    explicit qr(const pp::matrix& A);

    // Solve QRx = b  =>  Rx = Q^T b  (back-substitution)
    pp::vector solve(const pp::vector& b) const;

    // Determinant of the original matrix: det(Q) * product(R diagonal)
    double det() const;

    // Inverse of the original matrix
    pp::matrix inverse() const;
};

} // namespace pp
