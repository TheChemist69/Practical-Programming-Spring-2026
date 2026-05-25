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
    // Forward to the shared Linear_Equations vector implementation.
    #include "../Linear_Equations/vector.cc"
}
