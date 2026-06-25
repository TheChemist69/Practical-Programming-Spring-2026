// "bessel.h" header file.
// Bessel functions of integer order J_n(x) computed two ways:
//   (1) by complex contour integration around the origin (Part B1), and
//   (2) by the spectral trapezoid on the unit circle (Part C),
// plus a power-series reference for validation.
//
// Contour representation:
//   J_n(x) = 1/(2*pi*i) * oint_C z^(-n-1) exp( (x/2)(z - 1/z) ) dz ,
// where C encircles the origin once counter-clockwise.  Putting z = e^{i theta}
// turns this into the periodic real integral
//   J_n(x) = 1/(2*pi) * int_0^{2*pi} cos( n*theta - x*sin(theta) ) d(theta) ,
// whose integrand is periodic -> the trapezoid converges exponentially.
#pragma once

namespace pp {

// J_n(x) via contour integration around a regular polygon of `sides` vertices
// of the given `radius`, using the complex adaptive integrator.
double bessel_contour(int n, double x, int sides = 6, double radius = 1.0);

// J_n(x) via the periodic trapezoidal rule with N equally spaced samples
// (the spectral method on the unit circle).
double bessel_trapezoid(int n, double x, int N);

// Reference J_n(x) from the ascending power series (for validation).
double bessel_reference(int n, double x);

}  // namespace pp
