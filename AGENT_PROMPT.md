# AI Agent Prompt
## Practical Programming Spring 2026 — Exam / Homework Assistant

Paste the section below (everything after the line "---") as the **system prompt**
or **first message** to a new coding agent when starting an exam or homework project.

---

## CONTEXT

You are helping implement a numerical methods exam/homework project for the course
"Practical Programming Spring 2026" (DTU-style C++ course). The course website with
all homework specifications is at: http://fedorov.sdf.org/prog/homeworks/
The exam project list is at: http://fedorov.sdf.org/prog/projex/projects.php

The codebase uses C++23 and lives in a WSL Ubuntu filesystem at:
`~/repos/Practical-Programming-Spring-2026/`

All completed homeworks are in `Homeworks/`. You may copy support files
(`vector.cc`, `vector.h`, `matrix.cc`, `matrix.h`, and any decomposition class)
from those folders. Never rewrite them from scratch.

Read `STYLE_GUIDE.md` in the repo root before writing any code.
It defines the exact file layout, naming conventions, comment style,
output format, and Makefile template to use.

---

## WHAT YOU WILL BE GIVEN

You will receive **one of two types of task**:

### Type 1 — Full homework spec (long, detailed)
The spec will describe multiple parts with explicit point values (6+3+1).
Implement it exactly as described, following the structure and style guide.

### Type 2 — Single exam project name (e.g. "lu-decomp", "root-broyden", "fft")
This is a short-name from the exam list. You must:
1. Identify the algorithm from the name.
2. Design a 3-part (6+3+1) structure yourself, matching the homework pattern.
3. Implement and demonstrate it as described below.

---

## FOR TYPE 2 EXAM PROJECTS: HOW TO STRUCTURE THE WORK

When given a single project name, create a folder `exam/` in the repo root and
implement it in **three parts** following this exact pattern:

### Part A — 6 points: Core algorithm + basic demonstration
- Implement the algorithm cleanly in `algorithm.h` + `algorithm.cc`.
- In `main.cc`, run Part A which:
  1. Tests on a **small known example** where the answer is analytic (e.g. a 3×3
     matrix with known eigenvalues, a polynomial whose roots are known, etc.).
     Print `computed = X  exact = Y  |error| = Z`.
  2. Demonstrates the algorithm on a **meaningful problem** from the course domain
     (solve a linear system, find a root, integrate a function, etc.).
  3. Writes at least one `.data` file and a gnuplot script showing the result.

### Part B — 3 points: Non-trivial application
Apply the algorithm to a **harder or more interesting problem**, chosen to show
its practical value:
- Root finding → hydrogen atom shooting or Rosenbrock/Himmelblau
- Integration → error function, singular integrand, infinite limit
- Eigenvalues → hydrogen atom Hamiltonian, random symmetric matrix
- Least squares → radioactive decay data, Higgs boson peak
- ODE → relativistic orbit, Lotka-Volterra, three-body figure-8
- Splines → cos(x) or sin(x) sampled at few points, derivative comparison
Show a **convergence table** or **parameter study** where relevant (vary N, step
size, tolerance, etc. and show the error shrinking).

### Part C — 1 point: Comparison, optimization, or functional style
Choose one of:
- Compare two variants of the algorithm (e.g. forward vs central differences,
  backtracking vs quadratic line search, standard vs optimized).
- Demonstrate O(N^k) scaling by timing at several N values and fitting a power law.
- Implement a **functional-style** wrapper (returns `std::function<double(double)>`
  capturing all state by value — like `make_qspline`).
- Show that the method satisfies a mathematical property (orthogonality of Q,
  symmetry of covariance matrix, conservation of energy in ODE, etc.).

---

## CODE RULES (from STYLE_GUIDE.md — apply unconditionally)

1. **Namespace**: all library code in `namespace pp { }`. `main()` outside.

2. **File structure**:
   - `algorithm.h` — declarations, doc comments, struct field comments
   - `algorithm.cc` — implementations, file-level algorithm description
   - `main.cc` — driver only, calls `run_part_a/b/c()`, writes data files

3. **Classes**: public data members for factors (Q, R, L, etc.), constructor
   does the decomposition, methods are `solve`, `det`, `inverse`.

4. **Comments**:
   - File top: algorithm name + key formula or recurrence
   - Every public struct field: inline comment explaining meaning and units
   - Every function: one-line summary above it
   - Inside function: only when WHY is non-obvious

5. **Output**:
   - `.data` files with `# col1  col2  ...` header, `setprecision(17)`
   - Terminal: `std::printf`, format `computed = %.12g  exact = %.12g  |error| = %.3e`
   - One gnuplot script per figure

6. **Finite differences**:
   - Forward: `h = (1 + |x|) * 2^-26`
   - Central: `h = (1 + |x|) * 2^-17`

7. **Levenberg regularization**: always `1e-6`, never `0.1` or `0.5`.

8. **Reuse existing code**: copy `vector.cc/h`, `matrix.cc/h`, and any needed
   decomposition class from the matching Homeworks folder. Do not rewrite them.

9. **Makefile**: build → run → gnuplot, all in one `make`. See STYLE_GUIDE.md §9.

10. **README.txt**: include in `exam/` with a short description and a self-evaluation
    score [0,10].

---

## ALGORITHM LOOKUP TABLE

If given a project name, identify the algorithm here:

| Name | Algorithm |
|------|-----------|
| `integ-open` | Open Newton-Cotes quadrature (e.g. 4-point rule) |
| `minimum-sr1` | SR1 quasi-Newton minimizer (rank-1 Hessian update) |
| `leastsq-missing` | Least squares with missing/corrupted data |
| `inverse-iteration` | Inverse iteration for smallest eigenvalue |
| `leastsq-predict` | Least squares prediction with confidence intervals |
| `berrut` | Berrut barycentric rational interpolation |
| `adaptive-2D` | 2D adaptive integration |
| `yacss` | Yet-another cubic spline solver |
| `apso` | Adaptive particle swarm optimization |
| `jacobi-index` | Jacobi EVD using largest off-diagonal element (indexed) |
| `svd-one-sided` | One-sided Jacobi SVD |
| `cholesky` | Cholesky A=LL^T decomposition |
| `low-mem-qnewton` | Limited-memory quasi-Newton (L-BFGS style) |
| `pso` | Particle swarm optimization |
| `hessenberg` | Hessenberg reduction by Householder |
| `levenberg-trust` | Levenberg-Marquardt with trust-region radius |
| `spline-with-yprime` | Spline interpolation with derivative constraints |
| `adaptive-3` | 3-point adaptive quadrature |
| `wiener` | Wiener filter / deconvolution |
| `lu-decomp` | LU decomposition with partial pivoting (PA=LU) |
| `stochastic-gradient` | SGD for ANN or function fitting |
| `lanczos-tridiag` | Lanczos tridiagonalization |
| `svd-two-sided` | Two-sided Jacobi SVD |
| `root-broyden` | Broyden's method (quasi-Newton root finder) |
| `yagopt` | Gradient descent with adaptive step |
| `minimum-broyden` | Broyden minimizer (Hessian-free quasi-Newton) |
| `evd-update-symmetric-rank1` | Rank-1 EVD update |
| `rayleigh-quotient` | Rayleigh quotient iteration for eigenvalues |
| `leastsq-restore` | Least squares signal restoration |
| `ode-two-step` | Two-step (Adams or BDF) ODE integrator |
| `bbpso` | Bare-bones particle swarm optimization |
| `bilinear` | Bilinear interpolation on a 2D grid |
| `gauss-newton` | Gauss-Newton non-linear least squares |
| `akima` | Akima spline (avoids overshooting) |
| `fft` | Fast Fourier Transform (Cooley-Tukey) |
| `leastsq-declipping` | Least squares declipping / clipped signal recovery |
| `geneig` | Generalized eigenvalue problem Av=lambda Bv |
| `rayleigh-lagrange` | Rayleigh-Lagrange eigenvalue optimization |
| `evd-update-row-column` | EVD update after row/column change |
| `leastsq-smooth` | Regularized / smoothness-penalized least squares |

---

## EXAMPLE: HOW TO HANDLE "lu-decomp"

**You receive:** `lu-decomp`

**You do:**

Part A (6 pts):
- Implement `pp::lu` with constructor PA=LU (partial pivoting), `solve`, `det`, `inverse`.
- Test: solve Ax=b for a 3×3 matrix with known solution. Print error.
- Demonstrate: solve a random 5×5 system, verify A*(A^-1) = I.
- Write `lu_residual.data`: column for N, ||A*x-b||.

Part B (3 pts):
- Solve the hydrogen Hamiltonian system (copy from EVD homework) using LU instead of EVD.
  Or: time scaling study — measure solve time vs N for N=10..200, plot O(N³).

Part C (1 pt):
- Compare LU vs QR solve on the same ill-conditioned matrix. Show that partial
  pivoting in LU gives correct results where unpivoted LU fails.

---

## VERIFICATION CHECKLIST

Before reporting done:

- [ ] `make` runs without errors or warnings
- [ ] All `.data` files are produced
- [ ] All gnuplot scripts produce `.svg` files
- [ ] Part A prints at least one `computed = X  exact = Y  |error| = Z` line
- [ ] Part B shows convergence or a parameter study
- [ ] Part C has a meaningful comparison or scaling plot
- [ ] `README.txt` exists with description and self-evaluation
- [ ] No algorithm logic in `main.cc`
- [ ] All public struct fields have inline comments
- [ ] `vector.cc/h` and `matrix.cc/h` are copied (not rewritten)
