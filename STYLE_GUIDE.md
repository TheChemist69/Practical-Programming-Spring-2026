# Code Style Guide
## Practical Programming Spring 2026 — C++ Numerical Methods

This guide documents the conventions used across all homework projects.
Use it to ensure any new code (exam or otherwise) matches the established style.

---

## 1. Repository and File Layout

Every homework or exam project lives in its **own folder** under the repo root.

```
exam/
├── Makefile
├── README.txt          # short description + self-evaluation [0,10]
├── algorithm.h         # class/function declarations + doc comments
├── algorithm.cc        # implementations
├── main.cc             # driver: runs all parts, writes data files
├── plot_partA.gnuplot  # one gnuplot script per figure
└── plot_partB.gnuplot
```

Rules:
- One `.h` + one `.cc` per logical unit (never put implementation in the header).
- `main.cc` is **only** a driver — no algorithm logic there.
- Every output file has a corresponding gnuplot script that reads it.
- Compiled binaries (`main`, `*.o`) are produced by make and **not** committed.

---

## 2. Namespace

All library code lives in namespace `pp`:

```cpp
namespace pp {
    // your class/function here
} // namespace pp
```

`main()` itself is outside the namespace. The `run()` function that does the
actual work lives inside `pp` and is called from `main()`:

```cpp
namespace pp {
int run() {
    run_part_a();
    run_part_b();
    run_part_c();
    return 0;
}
} // namespace pp

int main() { return pp::run(); }
```

---

## 3. Vector and Matrix Classes

These classes are copied into each homework folder from the course examples.
**Never modify them.** Use them as-is.

```cpp
#include "vector.h"   // pp::vector — wraps std::vector<double>
#include "matrix.h"   // pp::matrix — stores columns as std::vector<pp::vector>
```

Key operations used throughout the homeworks:

```cpp
pp::vector v(n);          // zero-initialised length-n vector
v[i]                      // element access
v.norm()                  // Euclidean norm
v.dot(u)                  // dot product
v + u, v - u, v * scalar  // arithmetic (all return new vectors)

pp::matrix A(n, m);       // n×m zero matrix
A[j]                      // j-th column (a pp::vector&)
A[i, j]                   // element access (C++23 multi-index operator)
A.rows(), A.columns()
A.T()                     // transpose
pp::matrix::identity(n)   // n×n identity
```

When you need to write a new algorithm, model the interface after existing classes:

```cpp
class my_decomp {
public:
    pp::matrix Q;         // expose factors as public members (homework convention)
    pp::matrix R;

    explicit my_decomp(const pp::matrix& A); // constructor does the work

    pp::vector solve(const pp::vector& b) const;
    double det() const;
    pp::matrix inverse() const;
};
```

---

## 4. Header File Convention

```cpp
// "algorithm.h" header file.
// One-paragraph description of the algorithm and what it computes.
//
// Mathematical context:
//   A = Q*R  where Q has orthonormal columns, R is upper-triangular.
//
// Exposed interface:
//   solve(b)   -- solves Ax = b
//   det()      -- determinant
//   inverse()  -- matrix inverse
#pragma once

#include "matrix.h"   // always include what you need, nothing more

namespace pp {

// Structs: put field-level comments on every non-obvious field.
struct MyOptions {
    double acc = 1e-6;     // convergence threshold: stop when ||...|| < acc
    int max_iter = 100;    // hard limit on iterations
};

struct MyResult {
    pp::vector x;          // solution vector
    bool converged = false; // true if acc was reached
    int iterations = 0;    // how many steps were taken
};

// Function declarations: one-line summary above each.
// Solves Ax = b using the decomposition.
pp::vector my_solve(const pp::matrix& A, const pp::vector& b,
                    const MyOptions& opts = {});

} // namespace pp
```

---

## 5. Implementation File Convention

```cpp
// "algorithm.cc" implementation file.
// Algorithm name and brief mathematical description.
//
// Key formula or recurrence:
//   L[j,j] = sqrt(A[j,j] - sum_{k<j} L[j,k]^2)
//   L[i,j] = (A[i,j] - ...) / L[j,j]
#include "algorithm.h"
#include <cmath>
#include <stdexcept>

namespace pp {

// Constructor: one sentence on what it computes and any side effects.
my_decomp::my_decomp(const pp::matrix& A) {
    // implementation
}

// solve: one sentence. Omit obvious "this function solves..." prose.
pp::vector my_decomp::solve(const pp::vector& b) const {
    // implementation
}

} // namespace pp
```

Rules for comments in `.cc` files:
- Add a comment when the **why** is non-obvious (a subtle formula, a workaround, a numerical pitfall).
- Do **not** narrate what the code obviously does.
- Mathematical formulae are always worth a comment.
- One short line is better than a paragraph.

---

## 6. Main Driver Structure

Every project follows the **6 + 3 + 1 point** structure from the homework spec.
Structure `main.cc` accordingly even for exam projects:

```cpp
// "main.cc" — [Project name]: driver.
#include "algorithm.h"
#include <cstdio>
#include <cmath>
#include <fstream>

namespace pp {
namespace {

// Part A: core demonstration — show the algorithm works on simple cases.
void run_part_a() {
    std::puts("=== Part A (6 points): [core algorithm] ===\n");
    // 1. Debug/sanity test on a small known example.
    // 2. Main demonstration.
    // 3. Write data file(s) for plotting.
}

// Part B: extend to a harder or more interesting problem.
void run_part_b() {
    std::puts("=== Part B (3 points): [extension] ===\n");
}

// Part C: the 1-point task — comparison, optimization, or functional style.
void run_part_c() {
    std::puts("=== Part C (1 point): [extra] ===\n");
}

} // namespace

int run() {
    run_part_a();
    run_part_b();
    run_part_c();

    std::puts("Data files written:");
    std::puts("  partA.data");
    std::puts("  partB.data");
    return 0;
}

} // namespace pp

int main() { return pp::run(); }
```

---

## 7. Output and Data Files

- All numerical output goes to **plain text `.data` files**, one per plot.
- First line of every file is a gnuplot comment describing the columns:

```cpp
std::ofstream fout("result.data");
fout << "# x   f_numeric   f_exact   abs_error\n";
fout << std::setprecision(17);  // full double precision
for (...) {
    fout << x << ' ' << f_num << ' ' << f_exact << ' '
         << std::abs(f_num - f_exact) << '\n';
}
```

- Use `std::printf` for terminal output (not `std::cout`).
- Format terminal output so it is human-readable at a glance:

```cpp
std::printf("  computed = %.12g   exact = %.12g   |error| = %.3e\n",
            computed, exact, std::abs(computed - exact));
```

---

## 8. Gnuplot Scripts

One script per figure, named `plot_<what>.gnuplot`:

```gnuplot
set terminal svg size 800,500 background "white"
set output "result.svg"
set xlabel "x"
set ylabel "f(x)"
set title "Numerical vs exact"
set key top right
plot \
    "result.data" u 1:2 w l lw 2 t "numeric", \
    "result.data" u 1:3 w l lw 2 dt 2 t "exact"
```

The Makefile invokes gnuplot after building and running:

```makefile
all: main
    ./main
    gnuplot plot_partA.gnuplot
    gnuplot plot_partB.gnuplot
```

---

## 9. Makefile Template

```makefile
CC = g++
CFLAGS = -std=c++23 -O2 -Wall

SRC = main.cc algorithm.cc vector.cc matrix.cc
OBJ = $(SRC:.cc=.o)

main: $(OBJ)
    $(CC) $(CFLAGS) -o $@ $^

%.o: %.cc
    $(CC) $(CFLAGS) -c $@ $<

clean:
    rm -f *.o main *.data *.svg

.PHONY: clean
```

Copy `vector.cc`, `vector.h`, `matrix.cc`, `matrix.h` from the nearest homework
folder that uses the same interface — they are identical across all homeworks.

---

## 10. Error Handling

- Use `std::runtime_error` for fatal algorithmic failures (singular matrix, rank
  deficiency, non-positive definite).
- Return a status struct (converged, iterations, etc.) for soft failures
  (no convergence, max iterations reached) — never throw for these.
- Validate inputs only at the public interface boundary; trust internal calls.

```cpp
// Good: throw for a genuine precondition violation
if (A.rows() != A.columns())
    throw std::runtime_error("requires square matrix");

// Good: return status for iterative non-convergence
if (iter >= max_iter) {
    result.converged = false;
    return result;
}
```

---

## 11. Finite Differences

Use the correct step size for the method:

| Method | Optimal h | Exponent |
|--------|-----------|----------|
| Forward difference | `(1 + |x|) * 2^-26` | `-26` (≈ sqrt(eps)) |
| Central difference | `(1 + |x|) * 2^-17` | `-17` (≈ cbrt(eps)) |

```cpp
double h_forward = (1.0 + std::abs(x)) * std::pow(2.0, -26.0);
double h_central = (1.0 + std::abs(x)) * std::pow(2.0, -17.0);
```

Always use `max(|x|, 1)` as the scale so the step is non-zero when `x = 0`.

---

## 12. Numerical Method Patterns

### Newton root-finding
```
fx = f(x)
loop:
    build J by forward finite differences
    dx = QR(J).solve(-fx)
    line search: halve alpha until ||f(x+alpha*dx)|| < ||f(x)||
    x += alpha * dx
    stop if ||fx|| < acc or ||step|| < ||dx_fd||
```

### Newton minimization
```
loop:
    g = numerical gradient of phi
    H = numerical Hessian of phi
    dx = -(H + lm*I).solve(g)      // Levenberg: lm = 1e-6
    backtracking line search
    stop if ||g|| < acc
```

### Adaptive ODE driver (Dormand-Prince RK45)
```
loop:
    (yh, err) = rkstep45(f, x, y, h)
    tol = (acc + eps*||yh||) * sqrt(|h| / |b-a|)
    if err <= tol: accept step, advance x
    h *= min((tol/err)^0.2 * 0.95, 2.0)
```

### Adaptive quadrature
```
recursive(f, a, b, acc, eps):
    Q1 = low-order estimate on [a,b]
    Q2 = high-order estimate on [a,b]
    err = |Q1 - Q2|
    if err < acc + eps*|Q2|: return Q2
    return recursive(f, a, mid) + recursive(f, mid, b)
```

### QR solve pattern (used everywhere)
```cpp
pp::qr decomp(A);
pp::vector x = decomp.solve(b);
```

---

## 13. Reusable Classes Already in the Repo

When an exam project needs one of these, copy it from the matching homework:

| Need | Source folder |
|------|--------------|
| `pp::qr` (QR decomposition) | `Homeworks/Linear_Equations/` |
| `pp::lu` (LU decomposition) | `Homeworks/Linear_Equations/` |
| `pp::cholesky` | `Homeworks/Linear_Equations/` |
| `pp::EVD` (Jacobi eigenvalues) | `Homeworks/EVD/` |
| `rkstep45` + `driver` (ODE) | `Homeworks/OED/` |
| `pp::Ann` (neural network) | `Homeworks/Artificial Neural Networks/` |
| `pp::lsfit` (least squares) | `Homeworks/Least_Squares/` |
| `pp::qspline`, `pp::cspline` | `Homeworks/Splines/` |
| adaptive integrators | `Homeworks/Integration/` |
| Monte Carlo | `Homeworks/Monte Carlo integration/` |
| Newton solver | `Homeworks/Root_Finding/` |
| Newton minimizer | `Homeworks/Minimization/` |
