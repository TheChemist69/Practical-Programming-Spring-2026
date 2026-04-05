// "vector.cc" implementation file.
#include "vector.h"
#include <cmath>
#include <sstream>
#include <iostream>
#include <cstdio>

namespace pp {

// Construct a zero vector of given size
vector::vector(int n) : data(n) {}

// Construct from initializer list
vector::vector(std::initializer_list<NUMBER> list) : data(list) {}

// Size
int vector::size() const { return static_cast<int>(data.size()); }

// Resize
void vector::resize(int n) { data.resize(n); }

// Element access
NUMBER& vector::operator[](int i) { return data[i]; }
NUMBER  vector::operator[](int i) const { return data[i]; }

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
pp::vector vector::operator*(NUMBER s) const {
    pp::vector result(size());
    for (int i = 0; i < size(); i++) result[i] = data[i] * s;
    return result;
}

// Scalar division
pp::vector vector::operator/(NUMBER s) const {
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

pp::vector& vector::operator*=(NUMBER s) {
    for (int i = 0; i < size(); i++) data[i] *= s;
    return *this;
}

pp::vector& vector::operator/=(NUMBER s) {
    for (int i = 0; i < size(); i++) data[i] /= s;
    return *this;
}

// Dot product: u*v = sum_i( u_i v_i )
NUMBER vector::dot(const pp::vector& other) const {
    NUMBER sum = 0.0;
    for (int i = 0; i < size(); i++) sum += data[i] * other[i];
    return sum;
}

// Euclidean norm: ||v|| = sqrt(v*v)
NUMBER vector::norm() const {
    return std::sqrt(dot(*this));
}

// Print to stdout
void vector::print(std::string s) const {
    std::cout << s;
    for (int i = 0; i < size(); i++) printf("%9.3g ", (double)data[i]);
    printf("\n");
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
pp::vector operator*(NUMBER s, const pp::vector& v) {
    return v * s;
}

NUMBER dot(const pp::vector& a, const pp::vector& b) {
    return a.dot(b);
}

// Approximate equality: |a_i - b_i| <= acc + eps*max(|a_i|,|b_i|)  for all i
// Combines absolute tolerance (acc) and relative tolerance (eps).
bool approx(const pp::vector& a, const pp::vector& b, NUMBER acc, NUMBER eps) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); i++) {
        NUMBER diff = std::abs(a[i] - b[i]);
        NUMBER scale = std::max(std::abs(a[i]), std::abs(b[i]));
        if (diff > acc + eps * scale) return false;
    }
    return true;
}

// Map: apply function element-wise
pp::vector vector::map(std::function<NUMBER(NUMBER)> f) const {
    pp::vector result(size());
    for (int i = 0; i < size(); i++) result[i] = f(data[i]);
    return result;
}

} // namespace pp
