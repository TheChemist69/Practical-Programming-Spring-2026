// "vector.h" header file.
// Dense vector class reused across the root-finding homework.
#pragma once

#ifdef LONG_DOUBLE
	#define NUMBER long double
#else
	#define NUMBER double
#endif

#include <functional>
#include <initializer_list>
#include <string>
#include <vector>

namespace pp {

// Lightweight dynamic vector with element-wise arithmetic and norms.
class vector {
public:
    std::vector<NUMBER> data;

    vector() = default;
    explicit vector(int n);
    vector(std::initializer_list<NUMBER> list);

    int size() const;
    void resize(int n);

    NUMBER& operator[](int i);
    NUMBER operator[](int i) const;

    pp::vector operator-() const;
    pp::vector operator+(const pp::vector& other) const;
    pp::vector operator-(const pp::vector& other) const;
    pp::vector operator*(NUMBER s) const;
    pp::vector operator/(NUMBER s) const;

    pp::vector& operator+=(const pp::vector& other);
    pp::vector& operator-=(const pp::vector& other);
    pp::vector& operator*=(NUMBER s);
    pp::vector& operator/=(NUMBER s);

    NUMBER dot(const pp::vector& other) const;
    NUMBER norm() const;

    pp::vector map(std::function<NUMBER(NUMBER)> f) const;

    void print(std::string s = "") const;
    std::string to_string() const;
};

pp::vector operator*(NUMBER s, const pp::vector& v);
NUMBER dot(const pp::vector& a, const pp::vector& b);
bool approx(const pp::vector& a, const pp::vector& b, NUMBER acc = 1e-9, NUMBER eps = 1e-9);

} // namespace pp
