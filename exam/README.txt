================================================================
Adaptive Integration of Complex-Valued Functions
Exam project — Practical Programming, Spring 2026
Nikolaos Antonios Iakynthos Nemet
This project was created using Claude Code
================================================================

----------------------------------------------------------------
GOAL
----------------------------------------------------------------
Generalise an adaptive numerical integrator to integrate a complex-valued
function f(z) of a complex variable z along a straight line between two points
a and b in the complex plane, and to evaluate interesting contour integrals —
including the Bessel function J_n(x) of integer order via its contour-integral
representation.

The integrator is the course book's recursive open-4-point adaptive rule , carried
over from real to complex arithmetic.  It works throughout in
std::complex<double>: the error estimate |Q - q| is the complex modulus, a
straight contour is parametrised as z(t) = a + t(b - a), and a closed contour
is a polygon of such segments.  A real integral is simply one whose imaginary
part is zero, and its real part is taken at the end.

The work is organised in three parts (6 + 3 + 1 points).  Each part below
states its objective, how it is done, and the figures that show the result.


================================================================
PART A — A complex-valued adaptive contour integrator        (6 pts)
================================================================
Objective
  Build the engine: integration of a complex f(z) along a straight segment and
  around closed polygonal contours, with automatic accuracy control.

How it is done
  adaptive_unit() implements the open-4-point recursive rule on a real
  parameter t; integrate_segment() applies it to z(t) = a + t(b - a); and
  integrate_polygon() sums the segment integrals around a closed loop.  The
  engine is verified three ways:
    - it reproduces analytic antiderivatives,  ∫ f(z) dz = F(b) - F(a);
    - the result is path-independent for an analytic integrand;
    - it satisfies the residue theorem,  ∮ dz/z = 2πi   (and ∮ z² dz = 0).

Figures
  Part A has no plot: its correctness is shown numerically in the console
  report and re-checked by the built-in self-test (the residue theorem holds
  to ~1e-13).


================================================================
PART B — Real integrals evaluated by contour methods         (3 pts)
================================================================
Theme: the classic strength of contour integration — turning hard REAL
integrals into contour integrals that the engine evaluates directly.

B1 — Bessel functions J_n(x)
  Objective:  compute J_n(x) of integer order (the assignment's named example).
  How it is done:  J_n(x) = 1/(2πi) ∮ z^(-n-1) exp((x/2)(z - 1/z)) dz is
         evaluated around a polygon enclosing the origin, and validated against
         the power series and the recurrence
         J_{n-1}(x) + J_{n+1}(x) = (2n/x) J_n(x).
  Figures:
    plot_bessel.svg
      The contour points for J0, J1 and J2 sit exactly on the exact series
      curves, showing that the contour integral reproduces the Bessel functions
      correctly across the whole range of x.
    plot_bessel_convergence.svg
      The error of J0(2.5) stays at about 1e-17 for every polygon from a triangle
      to a 24-gon — the small wiggle is only rounding noise — so the result does
      not depend on the contour shape.  This is Cauchy's theorem: any polygon
      enclosing the origin gives the same value, which is why straight segments
      are enough to carry out a "circular" contour integral.

B2 — Trigonometric integrals over the unit circle
  Objective:  evaluate  ∫_0^{2π} dθ / (a + b cosθ) = 2π / √(a² - b²).
  How it is done:  the substitution z = e^{iθ}  (dθ = dz/(i z),
         cosθ = (z + 1/z)/2)  turns it into a contour integral on |z| = 1,
         evaluated on a 128-gon.  Pushing b toward a shows the adaptive rule
         holding accuracy while its evaluation count climbs as the pole nears
         the contour.
  Figures:
    plot_real_trig.svg
      The contour values land on the exact 2π/√(4 - b²) curve, including the steep
      rise as b → 2, confirming that the unit-circle contour integral really equals
      the real trigonometric integral and stays accurate even as the pole nears
      the contour.

B3 — A real-line integral by a rectangular contour
  Objective:  evaluate  ∫_{-∞}^{∞} dx / (1 + x²) = π.
  How it is done:  the contour is closed as a rectangle enclosing the pole
         z = i.  The residue theorem makes the whole loop equal 2πi·Res = π
         exactly; the bottom edge is the real integral itself, and it tends to
         π as the rectangle grows, because the other three edges vanish for
         large |z|.
  Figures:
    plot_real_rectangle.svg
      The gap |π - (bottom-edge integral)| falls steadily as the rectangle
      half-width R grows: the three non-real edges vanish as the rectangle
      grows, so the real-line integral converges to the exact contour value π.


================================================================
PART C — Spectral accuracy on a periodic contour             (1 pt)
================================================================
Objective
  Show when a simple rule beats the adaptive one.

How it is done
  On the circle z = e^{iθ} the Bessel integrand is periodic, and for a periodic
  integrand the plain trapezoidal rule converges exponentially.  J_3(5) is
  computed by the periodic trapezoid and by the book's (low-order) open-4
  adaptive rule on the same integral, and the error is compared against the
  number of integrand evaluations.  This is an honest "right tool for the job"
  comparison, not a claim that adaptivity is bad in general.

Figures
  plot_spectral.svg
      The trapezoid's error plunges to about 1e-16 by roughly 24 evaluations,
      while the adaptive rule's error falls slowly and needs far more evaluations
      to catch up.  For a smooth periodic integrand the trapezoid is spectrally
      accurate, so here it is both the right tool and far cheaper than general
      adaptivity.


================================================================
CONTENTS OF THE PROJECT
================================================================
cintegrator.h / cintegrator.cc — the adaptive integrator (complex)
    adaptive_unit(g, a, b)         adaptive integral of a complex g over [a, b]
    integrate_segment(f, a, b)     complex line integral along the segment a → b
    integrate_polygon(f, vertices) complex integral around a closed polygon
    integrate_rec(...)             (internal) the recursive open-4-point core
    CResult, COptions, IntegrationStatus   result / options / status types

bessel.h / bessel.cc — Bessel functions of integer order
    bessel_contour(n, x, sides, radius)  J_n(x) by contour integration (B1)
    bessel_trapezoid(n, x, N)            J_n(x) by the periodic trapezoid (C)
    bessel_reference(n, x)               J_n(x) by the power series (reference)

main.cc — orchestration only (no algorithm logic)
    run_part_a()        Part A: builds and verifies the integrator
    run_part_b()        Part B: Bessel, trig, and real-line integrals
    run_part_c()        Part C: spectral trapezoid vs the adaptive rule
    run_checks()        automated self-test of the headline results
    unit_circle(), trig_integrand(), one_over_1_plus_z2()   shared integrands
    emit() + LogScope   mirror all console output to Out.txt (RAII)

plot_*.gnuplot          one gnuplot script per figure (named in the parts above)
Out.txt                 the full console report, written on every run


================================================================
RESULTS
================================================================
  analytic antiderivatives, residue theorem : |err| ~ 1e-13
  Bessel J_n(x) vs series reference          : |err| ~ 1e-16
  polygon-shape independence (3- to 24-gon)  : agree to ~1e-17
  unit-circle trig integrals vs closed form  : |err| ~ 1e-12
  rectangle loop = π (residue theorem)       : |err| ~ 1e-11
  real-line integral: bottom edge → π        : edge error ~ 2/R
  spectral trapezoid J_3(5)                  : machine precision at N = 24
  automated self-test                        : 7/7 checks pass

All seven results above are re-checked by the built-in self-test on every run.


================================================================
SELF-EVALUATION
================================================================
Score: 9 / 10

The required task (complex adaptive contour integration with the Bessel
example) is complete, verified to machine precision, and cleanly structured
into separate compilation units with a single orchestrating main.  Beyond the
requirement it covers the classic uses of contour integration, trigonometric
integrals over the unit circle, a real-line integral by the residue theorem,
and (Part C) spectral accuracy on a periodic contour, stated honestly as a
"right tool for the job" comparison against the book's low-order rule rather
than a blanket claim.  Robustness is shown by a near-singular case (accuracy
held, cost rising), the error estimate is a conservative triangle-inequality
bound, and an automated self-test guards every headline result.  Held back from
10 only because a higher-order rule (Gauss-Kronrod) was not added; the project
deliberately generalises the book's open-4 rule instead, to stay grounded in
the course material and free of tabulated magic constants.
