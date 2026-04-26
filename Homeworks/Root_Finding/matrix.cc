// "matrix.cc" implementation file.
// Dense matrix storage and basic operations.

#include "matrix.h"

#include <cassert>

namespace pp {

matrix::matrix(int rows, int cols) : m_rows(rows), m_cols(cols), m_data(static_cast<std::size_t>(rows * cols), 0.0) {}

void matrix::resize(int rows, int cols) {
    m_rows = rows;
    m_cols = cols;
    m_data.assign(static_cast<std::size_t>(rows * cols), 0.0);
}

int matrix::rows() const { return m_rows; }
int matrix::cols() const { return m_cols; }

double& matrix::operator()(int i, int j) {
    assert(i >= 0 && i < m_rows);
    assert(j >= 0 && j < m_cols);
    return m_data[static_cast<std::size_t>(i * m_cols + j)];
}

double matrix::operator()(int i, int j) const {
    assert(i >= 0 && i < m_rows);
    assert(j >= 0 && j < m_cols);
    return m_data[static_cast<std::size_t>(i * m_cols + j)];
}

void matrix::set_zero() {
    for (double& x : m_data) x = 0.0;
}

vector matrix::column(int j) const {
    vector col(m_rows);
    for (int i = 0; i < m_rows; i++) {
        col[i] = (*this)(i, j);
    }
    return col;
}

void matrix::set_column(int j, const vector& col) {
    assert(col.size() == m_rows);
    for (int i = 0; i < m_rows; i++) {
        (*this)(i, j) = col[i];
    }
}

vector operator*(const matrix& A, const vector& x) {
    assert(A.cols() == x.size());
    vector y(A.rows());
    for (int i = 0; i < A.rows(); i++) {
        double sum = 0.0;
        for (int j = 0; j < A.cols(); j++) {
            sum += A(i, j) * x[j];
        }
        y[i] = sum;
    }
    return y;
}

matrix transpose(const matrix& A) {
    matrix T(A.cols(), A.rows());
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            T(j, i) = A(i, j);
        }
    }
    return T;
}

} // namespace pp
