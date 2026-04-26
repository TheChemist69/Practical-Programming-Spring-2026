// "qr.h" header file.
// QR decomposition and linear-system solve routines.
#pragma once

#include "matrix.h"
#include "vector.h"

namespace pp {

struct QRDecomposition {
    matrix Q;
    matrix R;
    bool success = true;
};

QRDecomposition qr_decompose(const matrix& A);
vector qr_solve(const QRDecomposition& qr, const vector& b, bool* success = nullptr);
vector solve_linear_system_qr(const matrix& A, const vector& b, bool* success = nullptr);

} // namespace pp
