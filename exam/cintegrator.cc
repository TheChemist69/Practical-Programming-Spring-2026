// "cintegrator.cc" implementation file.
// The recursive open-4-point adaptive engine (complex-valued) and the contour
// helpers built on it.
//
// Core rule on an interval [a,b] of a real parameter t (open: the endpoints are
// never evaluated, so integrable end singularities are tolerated):
//   nodes:  a + h/6, a + 2h/6, a + 4h/6, a + 5h/6        (h = b - a)
//   Q (higher order) = (2 f1 + f2 + f3 + 2 f4)/6 * h
//   q (lower order)  = (  f1 + f2 + f3 +   f4)/4 * h
//   error = |Q - q|  (the complex modulus)
//   accept if error <= acc + eps*|Q|, else split the interval, reuse (f1,f2) on
//   the left and (f3,f4) on the right, and rescale acc -> acc/sqrt(2).
#include "cintegrator.h"

#include <algorithm>
#include <cmath>
#include <optional>

namespace pp {
namespace {

// Node returned while descending the recursion tree.
struct Rec {
    cplx value{};
    std::size_t evaluations = 0;
    int depth = 0;
};

// Recursive core.  f2, f3 are the two inner-node values inherited from the
// parent call (empty on the first call for this interval).
Rec integrate_rec(const std::function<cplx(double)>& g,
                  double a, double b, double acc, double eps,
                  std::optional<cplx> f2, std::optional<cplx> f3,
                  int depth, const COptions& opt) {
    Rec out;
    out.depth = depth;

    const double h = b - a;
    const double x2 = a + 2.0 * h / 6.0;
    const double x3 = a + 4.0 * h / 6.0;

    // Reuse inherited inner nodes when present, else evaluate them.
    if (!f2) { f2 = g(x2); out.evaluations++; }
    if (!f3) { f3 = g(x3); out.evaluations++; }

    // The two outer open nodes are always fresh for this interval.
    const double x1 = a + h / 6.0;
    const double x4 = a + 5.0 * h / 6.0;
    const cplx f1 = g(x1);
    const cplx f4 = g(x4);
    out.evaluations += 2;

    // Embedded pair: higher-order Q and lower-order q.
    const cplx Q = (f1 * 2.0 + *f2 + *f3 + f4 * 2.0) * (h / 6.0);
    const cplx q = (f1 + *f2 + *f3 + f4) * (h / 4.0);
    const double err = std::abs(Q - q);
    const double tol = acc + eps * std::abs(Q);

    out.value = Q;

    // Stop on a non-finite value (also prevents infinite recursion on NaN),
    // on meeting the local tolerance, or on hitting the recursion guard.
    if (!std::isfinite(std::abs(Q)) || !std::isfinite(err)) return out;
    if (err <= tol) return out;
    if (depth >= opt.max_depth) return out;

    // Subdivide; distribute the absolute goal as acc/sqrt(2) per half.
    const double mid = 0.5 * (a + b);
    const double child_acc = acc / std::sqrt(2.0);

    // Left half inherits (f1,f2) as its inner nodes; right half inherits (f3,f4).
    const Rec left = integrate_rec(g, a, mid, child_acc, eps, f1, *f2, depth + 1, opt);
    const Rec right = integrate_rec(g, mid, b, child_acc, eps, *f3, f4, depth + 1, opt);

    out.value = left.value + right.value;
    out.evaluations += left.evaluations + right.evaluations;
    out.depth = std::max(left.depth, right.depth);
    return out;
}

}  // namespace

CResult adaptive_unit(const std::function<cplx(double)>& g,
                      double a, double b, const COptions& opt) {
    CResult r;
    if (a == b) return r;
    const Rec rec = integrate_rec(g, a, b, opt.acc, opt.eps, std::nullopt, std::nullopt, 0, opt);
    r.value = rec.value;
    r.evaluations = rec.evaluations;
    return r;
}

// integral_{a->b} f(z) dz, parametrized as z(t) = a + t(b-a), dz = (b-a) dt.
CResult integrate_segment(const std::function<cplx(cplx)>& f,
                          cplx a, cplx b, const COptions& opt) {
    const cplx d = b - a;
    const std::function<cplx(double)> g = [&](double t) -> cplx {
        return f(a + d * t) * d;
    };
    return adaptive_unit(g, 0.0, 1.0, opt);
}

// Closed polygon: sum the segment integrals over each edge, including the
// closing edge from the last vertex back to the first.
CResult integrate_polygon(const std::function<cplx(cplx)>& f,
                          const std::vector<cplx>& vertices,
                          const COptions& opt) {
    CResult total;
    const std::size_t n = vertices.size();
    if (n < 2) return total;

    for (std::size_t i = 0; i < n; i++) {
        const cplx a = vertices[i];
        const cplx b = vertices[(i + 1) % n];
        const CResult seg = integrate_segment(f, a, b, opt);
        total.value += seg.value;
        total.evaluations += seg.evaluations;
    }
    return total;
}

}  // namespace pp
