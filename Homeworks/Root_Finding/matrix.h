// "matrix.h" header file.
// Dense matrix utilities used for Jacobians and linear solves.
#pragma once

#include <vector>

#include "vector.h"

namespace pp {

class matrix {
private:
    int m_rows = 0;
    int m_cols = 0;
    std::vector<double> m_data;

public:
    matrix() = default;
    matrix(int rows, int cols);

    void resize(int rows, int cols);

    int rows() const;
    int cols() const;

    double& operator()(int i, int j);
    double operator()(int i, int j) const;

    void set_zero();

    vector column(int j) const;
    void set_column(int j, const vector& col);
};

vector operator*(const matrix& A, const vector& x);
matrix transpose(const matrix& A);

} // namespace pp
