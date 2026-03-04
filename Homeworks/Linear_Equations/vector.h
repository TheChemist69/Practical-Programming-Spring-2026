// "vector.h" header file.
// Column-vector class for linear algebra.
// Only doubles, no templates needed.
#pragma once

#include <vector>
#include <initializer_list>
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
    explicit vector(int n);

    // Construct from initializer list
    vector(std::initializer_list<double> list);

    // Size
    int size() const;

    // Element access
    double& operator[](int i);
    double  operator[](int i) const;

    // Unary negation
    pp::vector operator-() const;

    // Arithmetic operators
    pp::vector operator+(const pp::vector& other) const;
    pp::vector operator-(const pp::vector& other) const;
    pp::vector operator*(double s) const;
    pp::vector operator/(double s) const;

    // Compound assignment operators
    pp::vector& operator+=(const pp::vector& other);
    pp::vector& operator-=(const pp::vector& other);
    pp::vector& operator*=(double s);
    pp::vector& operator/=(double s);

    // Dot product
    double dot(const pp::vector& other) const;

    // Euclidean norm
    double norm() const;

    // String representation
    std::string to_string() const;
};

// Scalar * vector (free function)
pp::vector operator*(double s, const pp::vector& v);

// Dot product as free function
double dot(const pp::vector& a, const pp::vector& b);

} // namespace pp
