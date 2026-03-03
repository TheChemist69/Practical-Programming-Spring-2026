// "vector.h" header file.
// Column-vector class for linear algebra.
// Only doubles, no templates needed.
#pragma once

#include <vector>
#include <cmath>
#include <initializer_list>
#include <sstream>
#include <string>

namespace pp {

// ───────────────────────────────────────────────────────────────────
// pp::vector — column-vector backed by std::vector<double>
// ───────────────────────────────────────────────────────────────────
class vector {
public:
    std::vector<double> data;

    // Default constructor
    vector() = default;

    // Construct a zero vector of given size
    explicit vector(int n) : data(n) {}

    // Construct from initializer list
    vector(std::initializer_list<double> list) : data(list) {}

    // Copy and move are generated automatically (value semantics, like doubles)

    // Size
    int size() const { return static_cast<int>(data.size()); }

    // Element access — returns a reference so v[i] = x works
    double& operator[](int i) { return data[i]; }
    double  operator[](int i) const { return data[i]; }

    // Unary negation
    pp::vector operator-() const {
        pp::vector result(size());
        for (int i = 0; i < size(); i++) result[i] = -data[i];
        return result;
    }

    // Vector addition
    pp::vector operator+(const pp::vector& other) const {
        pp::vector result(size());
        for (int i = 0; i < size(); i++) result[i] = data[i] + other[i];
        return result;
    }

    // Vector subtraction
    pp::vector operator-(const pp::vector& other) const {
        pp::vector result(size());
        for (int i = 0; i < size(); i++) result[i] = data[i] - other[i];
        return result;
    }

    // Scalar multiplication (vector * scalar)
    pp::vector operator*(double s) const {
        pp::vector result(size());
        for (int i = 0; i < size(); i++) result[i] = data[i] * s;
        return result;
    }

    // Scalar division (vector / scalar)
    pp::vector operator/(double s) const {
        pp::vector result(size());
        for (int i = 0; i < size(); i++) result[i] = data[i] / s;
        return result;
    }

    // Compound assignment operators
    pp::vector& operator+=(const pp::vector& other) {
        for (int i = 0; i < size(); i++) data[i] += other[i];
        return *this;
    }

    pp::vector& operator-=(const pp::vector& other) {
        for (int i = 0; i < size(); i++) data[i] -= other[i];
        return *this;
    }

    pp::vector& operator*=(double s) {
        for (int i = 0; i < size(); i++) data[i] *= s;
        return *this;
    }

    pp::vector& operator/=(double s) {
        for (int i = 0; i < size(); i++) data[i] /= s;
        return *this;
    }

    // Dot product
    double dot(const pp::vector& other) const {
        double sum = 0.0;
        for (int i = 0; i < size(); i++) sum += data[i] * other[i];
        return sum;
    }

    // Euclidean norm
    double norm() const {
        return std::sqrt(dot(*this));
    }

    // String representation
    std::string to_string() const {
        std::ostringstream oss;
        oss << "{ ";
        for (int i = 0; i < size(); i++) {
            oss << data[i];
            if (i + 1 < size()) oss << ", ";
        }
        oss << " }";
        return oss.str();
    }
};

// Scalar * vector (free function)
inline pp::vector operator*(double s, const pp::vector& v) {
    return v * s;
}

// Dot product as free function
inline double dot(const pp::vector& a, const pp::vector& b) {
    return a.dot(b);
}

} // namespace pp
