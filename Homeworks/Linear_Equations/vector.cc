// "vector.cc" implementation file.
#include "vector.h"
#include <cmath>
#include <sstream>

namespace pp {

// Construct a zero vector of given size
vector::vector(int n) : data(n) {}

// Construct from initializer list
vector::vector(std::initializer_list<double> list) : data(list) {}

// Size
int vector::size() const { return static_cast<int>(data.size()); }

// Element access
double& vector::operator[](int i) { return data[i]; }
double  vector::operator[](int i) const { return data[i]; }

// Unary negation
pp::vector vector::operator-() const {
    pp::vector result(size());
    for (int i = 0; i < size(); i++) result[i] = -data[i];
    return result;
}

// Vector addition
pp::vector vector::operator+(const pp::vector& other) const {
    pp::vector result(size());
    for (int i = 0; i < size(); i++) result[i] = data[i] + other[i];
    return result;
}

// Vector subtraction
pp::vector vector::operator-(const pp::vector& other) const {
    pp::vector result(size());
    for (int i = 0; i < size(); i++) result[i] = data[i] - other[i];
    return result;
}

// Scalar multiplication
pp::vector vector::operator*(double s) const {
    pp::vector result(size());
    for (int i = 0; i < size(); i++) result[i] = data[i] * s;
    return result;
}

// Scalar division
pp::vector vector::operator/(double s) const {
    pp::vector result(size());
    for (int i = 0; i < size(); i++) result[i] = data[i] / s;
    return result;
}

// Compound assignment operators
pp::vector& vector::operator+=(const pp::vector& other) {
    for (int i = 0; i < size(); i++) data[i] += other[i];
    return *this;
}

pp::vector& vector::operator-=(const pp::vector& other) {
    for (int i = 0; i < size(); i++) data[i] -= other[i];
    return *this;
}

pp::vector& vector::operator*=(double s) {
    for (int i = 0; i < size(); i++) data[i] *= s;
    return *this;
}

pp::vector& vector::operator/=(double s) {
    for (int i = 0; i < size(); i++) data[i] /= s;
    return *this;
}

// Dot product
double vector::dot(const pp::vector& other) const {
    double sum = 0.0;
    for (int i = 0; i < size(); i++) sum += data[i] * other[i];
    return sum;
}

// Euclidean norm
double vector::norm() const {
    return std::sqrt(dot(*this));
}

// String representation
std::string vector::to_string() const {
    std::ostringstream oss;
    oss << "{ ";
    for (int i = 0; i < size(); i++) {
        oss << data[i];
        if (i + 1 < size()) oss << ", ";
    }
    oss << " }";
    return oss.str();
}

// Free functions
pp::vector operator*(double s, const pp::vector& v) {
    return v * s;
}

double dot(const pp::vector& a, const pp::vector& b) {
    return a.dot(b);
}

} // namespace pp
