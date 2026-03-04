// "matrix.h" header file.
// Matrix class for linear algebra.
// Only doubles, no templates needed.
#pragma once

#include <vector>
#include <string>
#include "vector.h"

namespace pp {

// ───────────────────────────────────────────────────────────────────
// pp::matrix — stores columns as std::vector<pp::vector>
// Column-major layout for data locality.
// ───────────────────────────────────────────────────────────────────
class matrix {
public:
    std::vector<pp::vector> cols;

    // Default constructor
    matrix() = default;

    // Construct an n×m zero matrix
    matrix(int n, int m);

    // Number of rows / columns
    int rows() const;
    int columns() const;

    // Unary negation
    pp::matrix operator-() const;

    // Column access: A[j] returns a reference to column j
    pp::vector&       operator[](int j);
    const pp::vector& operator[](int j) const;

    // Element access: A[i,j] (C++23 multidimensional subscript)
    double& operator[](int i, int j);
    double  operator[](int i, int j) const;

    // Static factory methods
    static pp::matrix identity(int n);
    static pp::matrix random(int n, int m);

    // Transpose
    pp::matrix transpose() const;

    // Matrix-vector multiplication
    pp::vector operator*(const pp::vector& v) const;

    // Matrix-matrix multiplication
    pp::matrix operator*(const pp::matrix& B) const;

    // Matrix addition / subtraction
    pp::matrix operator+(const pp::matrix& B) const;
    pp::matrix operator-(const pp::matrix& B) const;

    // Scalar multiplication / division
    pp::matrix operator*(double s) const;
    pp::matrix operator/(double s) const;

    // Compound assignment operators
    pp::matrix& operator+=(const pp::matrix& B);
    pp::matrix& operator-=(const pp::matrix& B);
    pp::matrix& operator*=(double s);
    pp::matrix& operator/=(double s);

    // Frobenius norm
    double norm() const;

    // String representation
    std::string to_string() const;
};

// Scalar * matrix (free function)
pp::matrix operator*(double s, const pp::matrix& A);

} // namespace pp
