
#include "vec.h"
#include <iostream>
#include <complex>

int main() {
	std::cout << "==== Basic Vector Operations (double) ====" << std::endl;
	vec<double> a(1.0, 2.0, 3.0);
	vec<double> b(4.0, 5.0, 6.0);
	vec<double> c = a + b;
	vec<double> d = a - b;
	vec<double> e = 2.0 * a;
	vec<double> f = b / 2.0;
	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl;
	std::cout << "c = a + b = " << c << std::endl;
	std::cout << "d = a - b = " << d << std::endl;
	std::cout << "e = 2.0 * a = " << e << std::endl;
	std::cout << "f = b / 2.0 = " << f << std::endl;

	std::cout << "\n==== Basic Vector Operations (float) ====" << std::endl;
	vec<float> af(1.0f, 2.0f, 3.0f);
	vec<float> bf(4.0f, 5.0f, 6.0f);
	std::cout << "af = " << af << std::endl;
	std::cout << "bf = " << bf << std::endl;
	std::cout << "af + bf = " << (af + bf) << std::endl;

	std::cout << "\n==== Basic Vector Operations (std::complex<double>) ====" << std::endl;
	using cd = std::complex<double>;
	vec<cd> acd(cd(1,1), cd(2,0), cd(0,2));
	vec<cd> bcd(cd(0,1), cd(1,1), cd(2,2));
	std::cout << "acd = " << acd << std::endl;
	std::cout << "bcd = " << bcd << std::endl;
	std::cout << "acd + bcd = " << (acd + bcd) << std::endl;

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
	vec<double> u(1, 0, 0), v(0, 1, 0), w(0, 0, 1);
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
	vec<double> a2 = a; // a2 is a copy of a
	std::cout << "Comparing a and a2 (copy of a):" << std::endl;
	std::cout << "a  = " << a << std::endl;
	std::cout << "a2 = " << a2 << std::endl;
	std::cout << "approx(a, a2): checks if a and a2 are approximately equal." << std::endl;
	std::cout << "Result: " << (approx(a, a2) ? "true" : "false") << std::endl;

	vec<double> a3 = a;
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
