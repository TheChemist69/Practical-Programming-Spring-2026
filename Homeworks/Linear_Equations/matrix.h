// "matrix.h" header file.
// Matrix class for linear algebra.
// Only doubles, no templates needed.
#pragma once

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <random>
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
    matrix(int n, int m) : cols(m, pp::vector(n)) {}

    // Copy and move are generated automatically (value semantics)

    // Number of rows
    int rows() const {
        if (cols.empty()) return 0;
        return cols[0].size();
    }

    // Number of columns
    int columns() const {
        return static_cast<int>(cols.size());
    }

    // Unary negation
    pp::matrix operator-() const {
        int n = rows();
        int m = columns();
        pp::matrix C(n, m);
        for (int j = 0; j < m; j++)
            for (int i = 0; i < n; i++)
                C[i, j] = -(*this)[i, j];
        return C;
    }

    // Column access: A[j] returns a reference to column j
    pp::vector&       operator[](int j)       { return cols[j]; }
    const pp::vector& operator[](int j) const { return cols[j]; }

    // Element access: A[i,j] (C++23 multidimensional subscript)
    double& operator[](int i, int j) { return cols[j][i]; }
    double  operator[](int i, int j) const { return cols[j][i]; }

    // Identity matrix
    static pp::matrix identity(int n) {
        pp::matrix I(n, n);
        for (int i = 0; i < n; i++) I[i, i] = 1.0;
        return I;
    }

    // Random matrix with entries in [0, 1)
    static pp::matrix random(int n, int m) {
        static std::mt19937 gen(42);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        pp::matrix A(n, m);
        for (int j = 0; j < m; j++)
            for (int i = 0; i < n; i++)
                A[i, j] = dist(gen);
        return A;
    }

    // Transpose
    pp::matrix transpose() const {
        int n = rows();
        int m = columns();
        pp::matrix T(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                T[i, j] = (*this)[j, i];
        return T;
    }

    // Matrix-vector multiplication: A * v
    pp::vector operator*(const pp::vector& v) const {
        int n = rows();
        int m = columns();
        pp::vector result(n);
        for (int j = 0; j < m; j++)
            for (int i = 0; i < n; i++)
                result[i] += (*this)[i, j] * v[j];
        return result;
    }

    // Matrix-matrix multiplication: A * B
    pp::matrix operator*(const pp::matrix& B) const {
        int n = rows();
        int m = B.columns();
        int p = columns();
        pp::matrix C(n, m);
        for (int j = 0; j < m; j++)
            for (int k = 0; k < p; k++)
                for (int i = 0; i < n; i++)
                    C[i, j] += (*this)[i, k] * B[k, j];
        return C;
    }

    // Matrix addition
    pp::matrix operator+(const pp::matrix& B) const {
        int n = rows();
        int m = columns();
        pp::matrix C(n, m);
        for (int j = 0; j < m; j++)
            for (int i = 0; i < n; i++)
                C[i, j] = (*this)[i, j] + B[i, j];
        return C;
    }

    // Matrix subtraction
    pp::matrix operator-(const pp::matrix& B) const {
        int n = rows();
        int m = columns();
        pp::matrix C(n, m);
        for (int j = 0; j < m; j++)
            for (int i = 0; i < n; i++)
                C[i, j] = (*this)[i, j] - B[i, j];
        return C;
    }

    // Scalar multiplication
    pp::matrix operator*(double s) const {
        int n = rows();
        int m = columns();
        pp::matrix C(n, m);
        for (int j = 0; j < m; j++)
            for (int i = 0; i < n; i++)
                C[i, j] = (*this)[i, j] * s;
        return C;
    }

    // Scalar division
    pp::matrix operator/(double s) const {
        int n = rows();
        int m = columns();
        pp::matrix C(n, m);
        for (int j = 0; j < m; j++)
            for (int i = 0; i < n; i++)
                C[i, j] = (*this)[i, j] / s;
        return C;
    }

    // Compound assignment operators
    pp::matrix& operator+=(const pp::matrix& B) {
        for (int j = 0; j < columns(); j++)
            for (int i = 0; i < rows(); i++)
                (*this)[i, j] += B[i, j];
        return *this;
    }

    pp::matrix& operator-=(const pp::matrix& B) {
        for (int j = 0; j < columns(); j++)
            for (int i = 0; i < rows(); i++)
                (*this)[i, j] -= B[i, j];
        return *this;
    }

    pp::matrix& operator*=(double s) {
        for (int j = 0; j < columns(); j++)
            for (int i = 0; i < rows(); i++)
                (*this)[i, j] *= s;
        return *this;
    }

    pp::matrix& operator/=(double s) {
        for (int j = 0; j < columns(); j++)
            for (int i = 0; i < rows(); i++)
                (*this)[i, j] /= s;
        return *this;
    }

    // Frobenius norm (useful for checking residuals)
    double norm() const {
        double sum = 0.0;
        for (int j = 0; j < columns(); j++)
            for (int i = 0; i < rows(); i++)
                sum += (*this)[i, j] * (*this)[i, j];
        return std::sqrt(sum);
    }

    // String representation
    std::string to_string() const {
        int n = rows();
        int m = columns();
        std::ostringstream oss;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                oss << (*this)[i, j];
                if (j + 1 < m) oss << "\t";
            }
            oss << "\n";
        }
        return oss.str();
    }
};

// Scalar * matrix (free function)
inline pp::matrix operator*(double s, const pp::matrix& A) {
    return A * s;
}

} // namespace pp
