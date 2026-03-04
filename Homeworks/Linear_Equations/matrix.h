// "matrix.h" header file.
// Matrix class for linear algebra.
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
    NUMBER  operator[](int i, int j) const;

    // Element access: A(i,j)
    NUMBER& operator()(int i, int j);
    NUMBER  operator()(int i, int j) const;

    // get/set element access
    NUMBER get(int i, int j) const;
    void set(int i, int j, NUMBER value);

    // Static factory methods
    static pp::matrix identity(int n);
    static pp::matrix random(int n, int m);

    // Transpose
    pp::matrix transpose() const;
    pp::matrix T() const;

    // Matrix-vector multiplication
    pp::vector operator*(const pp::vector& v) const;

    // Matrix-matrix multiplication
    pp::matrix operator*(const pp::matrix& B) const;

    // Matrix addition / subtraction
    pp::matrix operator+(const pp::matrix& B) const;
    pp::matrix operator-(const pp::matrix& B) const;

    // Scalar multiplication / division
    pp::matrix operator*(NUMBER s) const;
    pp::matrix operator/(NUMBER s) const;

    // Compound assignment operators
    pp::matrix& operator+=(const pp::matrix& B);
    pp::matrix& operator-=(const pp::matrix& B);
    pp::matrix& operator*=(NUMBER s);
    pp::matrix& operator/=(NUMBER s);

    // Frobenius norm
    NUMBER norm() const;

    // Print to stdout
    void print(std::string s = "") const;

    // String representation
    std::string to_string() const;
};

// Scalar * matrix (free function)
pp::matrix operator*(NUMBER s, const pp::matrix& A);

} // namespace pp
