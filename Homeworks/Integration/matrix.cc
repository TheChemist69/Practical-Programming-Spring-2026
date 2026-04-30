// "matrix.cc" implementation file.
#include "matrix.h"
// Forward to the shared Linear_Equations matrix implementation.
#include "../Linear_Equations/matrix.cc"
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

// Frobenius norm: ||A||_F = sqrt(sum_{ij}( A_{ij}^2 ))
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
