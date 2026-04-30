// "matrix.h" header file.
// Matrix class for linear algebra.
#pragma once

#include <vector>
#include <string>
#include "vector.h"

namespace pp {

// -----------------------------------------------------------------
// pp::matrix -- dense matrix for numerical linear algebra
//
// Stored in column-major order: each column is a pp::vector, and
// the matrix is a std::vector of those columns.  This layout gives
// good data locality for column-oriented algorithms (QR, etc.).
//
// Element access:
//   A[j]      -- reference to column j  (pp::vector&)
//   A[i,j]    -- element at row i, column j  (C++23 subscript)
//   A(i,j)    -- same, via operator()
//
// Supports matrix-vector (A*v) and matrix-matrix (A*B)
// multiplication, transpose (A.T()), Frobenius norm, and the
// usual element-wise arithmetic (+, -, scalar *, /).
//
// Factory methods:
//   matrix::identity(n) -- nxn identity matrix
//   matrix::random(n,m) -- nxm matrix with entries in [0,1)
// -----------------------------------------------------------------
class matrix {
public:
    std::vector<pp::vector> cols;

    // Default constructor
    matrix() = default;

    // Construct an nxm zero matrix
    matrix(int n, int m);

    // Number of rows / columns
    int rows() const;
    int columns() const;

    // Resize
    void resize(int n, int m);

    // Set to identity matrix (in-place)
    void setid();

    // Unary negation
    pp::matrix operator-() const;

    // Column access: A[j] returns a reference to column j
    pp::vector&       operator[](int j);
    const pp::vector& operator[](int j) const;

    // Element access: A[i,j] (C++23 multidimensional subscript)
    NUMBER& operator[](int i, int j);
    // Thin wrapper around the shared Linear_Equations matrix implementation.
    #pragma once

    #include "../Linear_Equations/matrix.h"
    NUMBER  operator()(int i, int j) const;
