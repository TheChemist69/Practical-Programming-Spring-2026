// "lu.h" header file.
// LU-decomposition via the Doolittle algorithm with partial pivoting.
// PA = LU where L is lower triangular (unit diagonal), U is upper triangular.
#pragma once

#include <vector>
#include "matrix.h"

namespace pp {

class lu {
public:
    pp::matrix L;
    pp::matrix U;
    std::vector<int> perm; // permutation vector
    int sign;              // sign of the permutation (+1 or -1)

    // Constructor: perform LU-decomposition of square matrix A
    explicit lu(const pp::matrix& A);

    // Forward substitution: solve Ly = b
    pp::vector forward_substitute(const pp::vector& b) const;

    // Back substitution: solve Ux = y
    pp::vector back_substitute(const pp::vector& y) const;

    // Solve Ax = b
    pp::vector solve(const pp::vector& b) const;

    // Determinant
    double det() const;

    // Inverse
    pp::matrix inverse() const;
};

} // namespace pp
