// "vector.h" header file.
// Column-vector class for linear algebra.
#pragma once

#ifdef LONG_DOUBLE
	#define NUMBER long double
#else
	#define NUMBER double
#endif

#include <vector>
#include <initializer_list>
#include <string>
#include <functional>

namespace pp {

// -----------------------------------------------------------------
// pp::vector -- dense column-vector for numerical linear algebra
//
// Backed by std::vector<NUMBER>, where NUMBER is double by default
// (or long double if LONG_DOUBLE is defined at compile time).
//
// Supports element-wise arithmetic (+, -, *, /), scalar
// multiplication from both sides, dot product, Euclidean norm
// ||v|| = sqrt(v*v), and element-wise function application via map().
//
// Element access uses 0-based indexing:  v[i].
// Construction from initializer list enables convenient literals:
//   pp::vector v = {1.0, 2.0, 3.0};
//
// Free functions provide commutativity (scalar * vector) and a
// standalone dot() as well as approximate-equality testing.
// -----------------------------------------------------------------
class vector {
public:
    std::vector<NUMBER> data;

    // Default constructor
    vector() = default;

    // Construct a zero vector of given size
    explicit vector(int n);

    // Construct from initializer list
    vector(std::initializer_list<NUMBER> list);

    // Size
    int size() const;

    // Resize
    void resize(int n);

    // Element access
    NUMBER& operator[](int i);
    NUMBER  operator[](int i) const;

    // Unary negation
    pp::vector operator-() const;

    // Arithmetic operators
    pp::vector operator+(const pp::vector& other) const;
    pp::vector operator-(const pp::vector& other) const;
    pp::vector operator*(NUMBER s) const;
    pp::vector operator/(NUMBER s) const;

    // Compound assignment operators
    pp::vector& operator+=(const pp::vector& other);
    pp::vector& operator-=(const pp::vector& other);
    pp::vector& operator*=(NUMBER s);
    pp::vector& operator/=(NUMBER s);

    // Dot product
    NUMBER dot(const pp::vector& other) const;

    // Euclidean norm
    NUMBER norm() const;

    // Apply a function element-wise
    pp::vector map(std::function<NUMBER(NUMBER)> f) const;

    // Print to stdout
    void print(std::string s = "") const;

    // String representation
    std::string to_string() const;
};

// Scalar * vector (free function)
pp::vector operator*(NUMBER s, const pp::vector& v);

// Dot product as free function
NUMBER dot(const pp::vector& a, const pp::vector& b);

// Approximate equality for vectors
bool approx(const pp::vector& a, const pp::vector& b, NUMBER acc = 1e-6, NUMBER eps = 1e-6);

} // namespace pp
