
#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <complex>
#include <type_traits>

template<typename T>
struct vec {
    T x, y, z;

    // constructors
    vec(T x, T y, T z) : x(x), y(y), z(z) {}
    vec() : vec(T{}, T{}, T{}) {}
    vec(const vec&) = default;
    vec(vec&&) = default;
    ~vec() = default;

    // assignment
    vec& operator=(const vec&) = default;
    vec& operator=(vec&&) = default;

    // arithmetic
    vec& operator+=(const vec& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
    vec& operator-=(const vec& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
    template<typename S>
    vec& operator*=(S s) { x *= s; y *= s; z *= s; return *this; }
    template<typename S>
    vec& operator/=(S s) { x /= s; y /= s; z /= s; return *this; }

    // utility
    void set(T a, T b, T c) { x = a; y = b; z = c; }
    void print(const std::string& s = "") const {
        if (!s.empty()) std::cout << s << ": ";
        std::cout << "[" << x << ", " << y << ", " << z << "]" << std::endl;
    }

    // vector math
    T dot(const vec& rhs) const { return x * rhs.x + y * rhs.y + z * rhs.z; }
    vec cross(const vec& rhs) const {
        return vec(
            y * rhs.z - z * rhs.y,
            z * rhs.x - x * rhs.z,
            x * rhs.y - y * rhs.x
        );
    }
    auto norm() const -> decltype(std::abs(x)) {
        using std::abs;
        return std::sqrt(abs(x) * abs(x) + abs(y) * abs(y) + abs(z) * abs(z));
    }

    // stream output
    friend std::ostream& operator<<(std::ostream& os, const vec& v) {
        os << "[" << v.x << ", " << v.y << ", " << v.z << "]";
        return os;
    }
};

// non-member operators
template<typename T>
vec<T> operator-(const vec<T>& v) { return vec<T>(-v.x, -v.y, -v.z); }

template<typename T>
vec<T> operator-(const vec<T>& a, const vec<T>& b) { return vec<T>(a.x - b.x, a.y - b.y, a.z - b.z); }

template<typename T>
vec<T> operator+(const vec<T>& a, const vec<T>& b) { return vec<T>(a.x + b.x, a.y + b.y, a.z + b.z); }

template<typename T, typename S>
vec<T> operator*(const vec<T>& v, S s) { return vec<T>(v.x * s, v.y * s, v.z * s); }

template<typename T, typename S>
vec<T> operator*(S s, const vec<T>& v) { return v * s; }

template<typename T, typename S>
vec<T> operator/(const vec<T>& v, S s) { return vec<T>(v.x / s, v.y / s, v.z / s); }

// approximate equality
template<typename T>
bool approx(const vec<T>& a, const vec<T>& b, double acc = 1e-6, double eps = 1e-6) {
    auto approx_scalar = [](auto x, auto y, double acc, double eps) {
        using std::abs;
        return abs(x - y) <= acc + eps * abs(x);
    };
    return approx_scalar(a.x, b.x, acc, eps)
        && approx_scalar(a.y, b.y, acc, eps)
        && approx_scalar(a.z, b.z, acc, eps);
}
