// "vector.cc" implementation file.
// Dense vector operations for numerical methods.

#include "vector.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>

namespace pp {

vector::vector(int n) : data(n) {}
vector::vector(std::initializer_list<NUMBER> list) : data(list) {}

int vector::size() const { return static_cast<int>(data.size()); }
void vector::resize(int n) { data.resize(n); }

NUMBER& vector::operator[](int i) { return data[i]; }
NUMBER vector::operator[](int i) const { return data[i]; }

pp::vector vector::operator-() const {
    pp::vector result(size());
    for (int i = 0; i < size(); i++) result[i] = -data[i];
    return result;
}

pp::vector vector::operator+(const pp::vector& other) const {
    pp::vector result(size());
    for (int i = 0; i < size(); i++) result[i] = data[i] + other[i];
    return result;
}

pp::vector vector::operator-(const pp::vector& other) const {
    pp::vector result(size());
    for (int i = 0; i < size(); i++) result[i] = data[i] - other[i];
    return result;
}

pp::vector vector::operator*(NUMBER s) const {
    pp::vector result(size());
    for (int i = 0; i < size(); i++) result[i] = data[i] * s;
    return result;
}

pp::vector vector::operator/(NUMBER s) const {
    pp::vector result(size());
    for (int i = 0; i < size(); i++) result[i] = data[i] / s;
    return result;
}

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

NUMBER vector::dot(const pp::vector& other) const {
    NUMBER sum = 0.0;
    for (int i = 0; i < size(); i++) sum += data[i] * other[i];
    return sum;
}

NUMBER vector::norm() const { return std::sqrt(dot(*this)); }

pp::vector vector::map(std::function<NUMBER(NUMBER)> f) const {
    pp::vector result(size());
    for (int i = 0; i < size(); i++) result[i] = f(data[i]);
    return result;
}

void vector::print(std::string s) const {
    std::cout << s;
    for (int i = 0; i < size(); i++) std::printf("%12.6g ", static_cast<double>(data[i]));
    std::printf("\n");
}

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

pp::vector operator*(NUMBER s, const pp::vector& v) { return v * s; }
NUMBER dot(const pp::vector& a, const pp::vector& b) { return a.dot(b); }

bool approx(const pp::vector& a, const pp::vector& b, NUMBER acc, NUMBER eps) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); i++) {
        const NUMBER diff = std::abs(a[i] - b[i]);
        const NUMBER scale = std::max(std::abs(a[i]), std::abs(b[i]));
        if (diff > acc + eps * scale) return false;
    }
    return true;
}

} // namespace pp
