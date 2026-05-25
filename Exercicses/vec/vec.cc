#include "vec.h"
#include <iostream>
#include <cmath>

// --- Member functions ---

double vec::dot(const vec& rhs) const {
	return x * rhs.x + y * rhs.y + z * rhs.z;
}

vec vec::cross(const vec& rhs) const {
	return vec(
		y * rhs.z - z * rhs.y,
		z * rhs.x - x * rhs.z,
		x * rhs.y - y * rhs.x
	);
}

double vec::norm() const {
	return std::sqrt(x * x + y * y + z * z);
}

vec& vec::operator+=(const vec& rhs) {
	x += rhs.x;
	y += rhs.y;
	z += rhs.z;
	return *this;
}

vec& vec::operator-=(const vec& rhs) {
	x -= rhs.x;
	y -= rhs.y;
	z -= rhs.z;
	return *this;
}

vec& vec::operator*=(double s) {
	x *= s;
	y *= s;
	z *= s;
	return *this;
}

vec& vec::operator/=(double s) {
	x /= s;
	y /= s;
	z /= s;
	return *this;
}

void vec::print(const std::string& s) const {
	if (!s.empty()) std::cout << s << ": ";
	std::cout << "[" << x << ", " << y << ", " << z << "]" << std::endl;
}

// --- Non-member functions ---

std::ostream& operator<<(std::ostream& os, const vec& v) {
	os << "[" << v.x << ", " << v.y << ", " << v.z << "]";
	return os;
}

vec operator-(const vec& v) {
	return vec(-v.x, -v.y, -v.z);
}

vec operator-(const vec& a, const vec& b) {
	return vec(a.x - b.x, a.y - b.y, a.z - b.z);
}

vec operator+(const vec& a, const vec& b) {
	return vec(a.x + b.x, a.y + b.y, a.z + b.z);
}

vec operator*(const vec& v, double s) {
	return vec(v.x * s, v.y * s, v.z * s);
}

vec operator*(double s, const vec& v) {
	return v * s;
}

vec operator/(const vec& v, double s) {
	return vec(v.x / s, v.y / s, v.z / s);
}

// --- Approximate equality ---

namespace {
	bool approx_double(double a, double b, double acc = 1e-6, double eps = 1e-6) {
		return std::fabs(a - b) <= acc + eps * std::fabs(a);
	}
}

bool approx(const vec& a, const vec& b, double acc, double eps) {
	return approx_double(a.x, b.x, acc, eps)
		&& approx_double(a.y, b.y, acc, eps)
		&& approx_double(a.z, b.z, acc, eps);
}
