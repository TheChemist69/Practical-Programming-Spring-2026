// "bessel.cc" implementation file.
// Integer-order Bessel functions via contour integration, via the spectral
// trapezoid, and via the power series (reference).
#include "bessel.h"

#include <cmath>
#include <numbers>
#include <vector>

#include "cintegrator.h"

namespace pp {

using std::numbers::pi;
using namespace std::complex_literals;  // the imaginary unit 1i

// J_n(x) = 1/(2 pi i) oint z^(-n-1) exp((x/2)(z - 1/z)) dz.
// The contour is a regular polygon (radius `radius`, `sides` edges) around the
// origin; by Cauchy's theorem the result is independent of the polygon shape.
double bessel_contour(int n, double x, int sides, double radius) {
    const auto integrand = [n, x](cplx z) -> cplx {
        return std::pow(z, -n - 1) * std::exp(0.5 * x * (z - 1.0 / z));
    };

    std::vector<cplx> vertices;
    vertices.reserve(sides);
    for (int k = 0; k < sides; k++) {
        const double theta = 2.0 * pi * k / sides;
        vertices.push_back(radius * cplx(std::cos(theta), std::sin(theta)));
    }

    COptions opt;
    opt.acc = 1e-11;
    opt.eps = 1e-11;
    const CResult r = integrate_polygon(integrand, vertices, opt);

    // Divide by 2*pi*i to recover J_n; the result is real up to rounding.
    const cplx J = r.value / (2.0 * pi * 1i);
    return J.real();
}

// J_n(x) = 1/(2 pi) int_0^{2 pi} cos(n theta - x sin theta) dtheta.
// For a periodic integrand the trapezoid is just the sample mean.
double bessel_trapezoid(int n, double x, int N) {
    double sum = 0.0;
    for (int k = 0; k < N; k++) {
        const double theta = 2.0 * pi * k / N;
        sum += std::cos(n * theta - x * std::sin(theta));
    }
    return sum / N;
}

// Ascending series: J_n(x) = sum_{m>=0} (-1)^m / (m! (m+n)!) * (x/2)^(2m+n).
// Accurate for the small/moderate x used here (x <= 12); at large x the
// alternating terms cancel and a few digits would be lost.
double bessel_reference(int n, double x) {
    const double half = 0.5 * x;

    // First term (m=0): (x/2)^n / n!
    double term = 1.0;
    for (int i = 1; i <= n; i++) term *= half / i;

    double sum = term;
    for (int m = 1; m < 100; m++) {
        // ratio term_m / term_{m-1} = -(x/2)^2 / (m (m+n))
        term *= -(half * half) / (static_cast<double>(m) * (m + n));
        sum += term;
        if (std::abs(term) < 1e-18 * std::abs(sum)) break;
    }
    return sum;
}

}  // namespace pp
