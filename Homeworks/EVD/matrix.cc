// "matrix.cc" implementation file.
#include "matrix.h"
#include <cmath>
#include <sstream>
#include <random>
#include <iostream>
#include <cstdio>

namespace pp {

// Construct an n×m zero matrix
matrix::matrix(int n, int m) : cols(m, pp::vector(n)) {}

// Number of rows
int matrix::rows() const {
    if (cols.empty()) return 0;
    return cols[0].size();
}

// Number of columns
int matrix::columns() const {
    return static_cast<int>(cols.size());
}

// Resize
void matrix::resize(int n, int m) {
    cols.resize(m, pp::vector(n));
    for (int j = 0; j < m; j++) cols[j].resize(n);
}

// Set to identity (in-place)
void matrix::setid() {
    int n = rows();
    int m = columns();
    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++)
            cols[j][i] = (i == j) ? 1.0 : 0.0;
}

// Unary negation
pp::matrix matrix::operator-() const {
    int n = rows();
    int m = columns();
    pp::matrix C(n, m);
    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++)
            C[i, j] = -(*this)[i, j];
    return C;
}

// Column access
pp::vector& matrix::operator[](int j) { return cols[j]; }
const pp::vector& matrix::operator[](int j) const { return cols[j]; }

// Element access (C++23)
NUMBER& matrix::operator[](int i, int j) { return cols[j][i]; }
NUMBER  matrix::operator[](int i, int j) const { return cols[j][i]; }

// Element access operator()
NUMBER& matrix::operator()(int i, int j) { return cols[j][i]; }
NUMBER  matrix::operator()(int i, int j) const { return cols[j][i]; }

// get/set
NUMBER matrix::get(int i, int j) const { return cols[j][i]; }
void matrix::set(int i, int j, NUMBER value) { cols[j][i] = value; }

// Identity matrix
pp::matrix matrix::identity(int n) {
    pp::matrix I(n, n);
    for (int i = 0; i < n; i++) I[i, i] = 1.0;
    return I;
}

// Random matrix with entries in [0, 1)
pp::matrix matrix::random(int n, int m) {
    static std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    pp::matrix A(n, m);
    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++)
            A[i, j] = dist(gen);
    return A;
}

// Transpose
pp::matrix matrix::transpose() const {
    int n = rows();
    int m = columns();
    pp::matrix T(m, n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            T[i, j] = (*this)[j, i];
    return T;
}

pp::matrix matrix::T() const { return transpose(); }

// Matrix-vector multiplication
pp::vector matrix::operator*(const pp::vector& v) const {
    int n = rows();
    int m = columns();
    pp::vector result(n);
    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++)
            result[i] += (*this)[i, j] * v[j];
    return result;
}

// Matrix-matrix multiplication
pp::matrix matrix::operator*(const pp::matrix& B) const {
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
pp::matrix matrix::operator+(const pp::matrix& B) const {
    int n = rows();
    int m = columns();
    pp::matrix C(n, m);
    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++)
            C[i, j] = (*this)[i, j] + B[i, j];
    return C;
}

// Matrix subtraction
pp::matrix matrix::operator-(const pp::matrix& B) const {
    int n = rows();
    int m = columns();
    pp::matrix C(n, m);
    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++)
            C[i, j] = (*this)[i, j] - B[i, j];
    return C;
}

// Scalar multiplication
pp::matrix matrix::operator*(NUMBER s) const {
    int n = rows();
    int m = columns();
    pp::matrix C(n, m);
    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++)
            C[i, j] = (*this)[i, j] * s;
    return C;
}

// Scalar division
pp::matrix matrix::operator/(NUMBER s) const {
    int n = rows();
    int m = columns();
    pp::matrix C(n, m);
    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++)
            C[i, j] = (*this)[i, j] / s;
    return C;
}

// Compound assignment operators
pp::matrix& matrix::operator+=(const pp::matrix& B) {
    for (int j = 0; j < columns(); j++)
        for (int i = 0; i < rows(); i++)
            (*this)[i, j] += B[i, j];
    return *this;
}

pp::matrix& matrix::operator-=(const pp::matrix& B) {
    for (int j = 0; j < columns(); j++)
        for (int i = 0; i < rows(); i++)
            (*this)[i, j] -= B[i, j];
    return *this;
}

pp::matrix& matrix::operator*=(NUMBER s) {
    for (int j = 0; j < columns(); j++)
        for (int i = 0; i < rows(); i++)
            (*this)[i, j] *= s;
    return *this;
}

pp::matrix& matrix::operator/=(NUMBER s) {
    for (int j = 0; j < columns(); j++)
        for (int i = 0; i < rows(); i++)
            (*this)[i, j] /= s;
    return *this;
}

// Frobenius norm
NUMBER matrix::norm() const {
    NUMBER sum = 0.0;
    for (int j = 0; j < columns(); j++)
        for (int i = 0; i < rows(); i++)
            sum += (*this)[i, j] * (*this)[i, j];
    return std::sqrt(sum);
}

// Print to stdout
void matrix::print(std::string s) const {
    std::cout << s << std::endl;
    for (int i = 0; i < rows(); i++) {
        for (int j = 0; j < columns(); j++)
            printf("%9.3g ", (double)(*this)[i, j]);
        printf("\n");
    }
    printf("\n");
}

// String representation
std::string matrix::to_string() const {
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

// Free function
pp::matrix operator*(NUMBER s, const pp::matrix& A) {
    return A * s;
}

} // namespace pp
