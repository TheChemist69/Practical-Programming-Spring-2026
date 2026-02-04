#include "vec.h"
#include <iostream>

int main() {
	std::cout << "==== Basic Vector Operations ====" << std::endl;
	vec a(1.0, 2.0, 3.0);
	vec b(4.0, 5.0, 6.0);
	vec c = a + b;
	vec d = a - b;
	vec e = 2.0 * a;
	vec f = b / 2.0;
	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl;
	std::cout << "c = a + b = " << c << std::endl;
	std::cout << "d = a - b = " << d << std::endl;
	std::cout << "e = 2.0 * a = " << e << std::endl;
	std::cout << "f = b / 2.0 = " << f << std::endl;

	std::cout << "\n==== Compound Assignment Demonstration ====" << std::endl;
	a += b;
	std::cout << "a += b: " << a << std::endl;
	b -= a;
	std::cout << "b -= a: " << b << std::endl;
	c *= 3.0;
	std::cout << "c *= 3.0: " << c << std::endl;
	d /= 2.0;
	std::cout << "d /= 2.0: " << d << std::endl;

	std::cout << "\n==== Print Method Demo ====" << std::endl;
	a.print("Debug a");

	std::cout << "\n==== Dot Product ====" << std::endl;
	vec u(1, 0, 0), v(0, 1, 0), w(0, 0, 1);
	std::cout << "u = " << u << ", v = " << v << ", w = " << w << std::endl;
	std::cout << "u . v = " << u.dot(v) << " (should be 0)" << std::endl;
	std::cout << "u . u = " << u.dot(u) << " (should be 1)" << std::endl;
	std::cout << "a . b = " << a.dot(b) << std::endl;

	std::cout << "\n==== Cross Product ====" << std::endl;
	std::cout << "u x v = " << u.cross(v) << " (should be [0, 0, 1])" << std::endl;
	std::cout << "v x w = " << v.cross(w) << " (should be [1, 0, 0])" << std::endl;
	std::cout << "a x b = " << a.cross(b) << std::endl;

	std::cout << "\n==== Norm (Magnitude) ====" << std::endl;
	std::cout << "|u| = " << u.norm() << " (should be 1)" << std::endl;
	std::cout << "|a| = " << a.norm() << std::endl;

	std::cout << "\n==== Approximate Equality ====" << std::endl;
	vec a2 = a; // a2 is a copy of a
	std::cout << "Comparing a and a2 (copy of a):" << std::endl;
	std::cout << "a  = " << a << std::endl;
	std::cout << "a2 = " << a2 << std::endl;
	std::cout << "approx(a, a2): checks if a and a2 are approximately equal." << std::endl;
	std::cout << "Result: " << (approx(a, a2) ? "true" : "false") << std::endl;

	vec a3 = a;
	a3.x += 1e-7;
	std::cout << "\nComparing a and a3 (a3.x += 1e-7):" << std::endl;
	std::cout << "a  = " << a << std::endl;
	std::cout << "a3 = " << a3 << std::endl;
	std::cout << "approx(a, a3): " << (approx(a, a3) ? "true" : "false") << std::endl;

	a3.x += 1e-3;
	std::cout << "\nComparing a and a3 (a3.x += 1e-3):" << std::endl;
	std::cout << "a  = " << a << std::endl;
	std::cout << "a3 = " << a3 << std::endl;
	std::cout << "approx(a, a3): " << (approx(a, a3) ? "true" : "false") << std::endl;

	return 0;
}
