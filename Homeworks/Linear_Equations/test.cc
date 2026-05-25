// test.cc — Unified test suite for the pp:: linear algebra library.
// Tests every operation and function with clear, readable output.
// Build:  g++ -std=c++23 -O2 -c test.cc -o test.o
// Link:   g++ test.o vector.o matrix.o qr.o lu.o cholesky.o householder.o givens.o -o test
// Run:    ./test
#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include "vector.h"
#include "matrix.h"
#include "qr.h"
#include "lu.h"
#include "cholesky.h"
#include "householder.h"
#include "givens.h"

// ── Globals ──────────────────────────────────────────────────────────
static const NUMBER tol = 1e-9;
static int total_pass = 0, total_fail = 0;

// ── Helpers ──────────────────────────────────────────────────────────
void check(bool ok, const char* desc) {
    if (ok) { std::cout << "  [PASS] " << desc << "\n"; total_pass++; }
    else    { std::cout << "  [FAIL] " << desc << "\n"; total_fail++; }
}

bool mat_approx(const pp::matrix& A, const pp::matrix& B) {
    return (A - B).norm() < tol;
}
bool vec_approx(const pp::vector& a, const pp::vector& b) {
    return (a - b).norm() < tol;
}
bool is_upper_tri(const pp::matrix& R) {
    for (int j = 0; j < R.columns(); j++)
        for (int i = j + 1; i < R.rows(); i++)
            if (std::abs(R[i, j]) > tol) return false;
    return true;
}
bool is_lower_tri(const pp::matrix& L) {
    for (int j = 0; j < L.columns(); j++)
        for (int i = 0; i < j; i++)
            if (std::abs(L[i, j]) > tol) return false;
    return true;
}
bool is_identity(const pp::matrix& M) {
    for (int i = 0; i < M.rows(); i++)
        for (int j = 0; j < M.columns(); j++)
            if (std::abs(M[i, j] - (i == j ? 1.0 : 0.0)) > tol) return false;
    return true;
}

void section(const char* title) {
    std::cout << "\n══════════════════════════════════════════════════════════\n";
    std::cout << "  " << title << "\n";
    std::cout << "══════════════════════════════════════════════════════════\n\n";
}

void subsection(const char* title) {
    std::cout << "  ── " << title << " ──\n";
}

// ═════════════════════════════════════════════════════════════════════
//  1. pp::vector — construction, element access, resize
// ═════════════════════════════════════════════════════════════════════
void test_vector_basics() {
    section("pp::vector — basics");

    // Construction
    subsection("Construction");
    pp::vector v0;
    check(v0.size() == 0, "Default constructor creates empty vector");

    pp::vector v3(4);
    check(v3.size() == 4, "vector(4) has size 4");
    check(v3[0] == 0 && v3[3] == 0, "vector(4) is zero-initialized");

    pp::vector vi = {1.0, 2.0, 3.0};
    check(vi.size() == 3, "Initializer list {1,2,3} has size 3");
    check(vi[0] == 1 && vi[1] == 2 && vi[2] == 3, "Initializer list values correct");

    // Element access
    subsection("Element access (operator[])");
    vi[1] = 42.0;
    check(vi[1] == 42.0, "Write vi[1]=42, read vi[1]==42");

    // Resize
    subsection("Resize");
    pp::vector vr = {10, 20, 30};
    vr.resize(5);
    check(vr.size() == 5, "resize(5) → size 5");
    check(vr[0] == 10 && vr[2] == 30, "resize preserves existing data");
    check(vr[3] == 0 && vr[4] == 0, "resize zero-fills new elements");
    vr.resize(2);
    check(vr.size() == 2, "resize(2) shrinks to size 2");
    check(vr[0] == 10 && vr[1] == 20, "shrink preserves surviving data");
}

// ═════════════════════════════════════════════════════════════════════
//  2. pp::vector — arithmetic operations
// ═════════════════════════════════════════════════════════════════════
void test_vector_arithmetic() {
    section("pp::vector — arithmetic");

    pp::vector a = {1, 2, 3};
    pp::vector b = {4, 5, 6};

    // Addition
    subsection("Addition (a + b)");
    pp::vector c = a + b;
    std::cout << "    a = " << a.to_string() << "\n";
    std::cout << "    b = " << b.to_string() << "\n";
    std::cout << "  a+b = " << c.to_string() << "\n";
    check(c[0] == 5 && c[1] == 7 && c[2] == 9, "a+b = {5,7,9}");

    // Subtraction
    subsection("Subtraction (a - b)");
    pp::vector d = a - b;
    std::cout << "  a-b = " << d.to_string() << "\n";
    check(d[0] == -3 && d[1] == -3 && d[2] == -3, "a-b = {-3,-3,-3}");

    // Unary negation
    subsection("Unary negation (-a)");
    pp::vector neg = -a;
    std::cout << "   -a = " << neg.to_string() << "\n";
    check(neg[0] == -1 && neg[1] == -2 && neg[2] == -3, "-a = {-1,-2,-3}");

    // Scalar multiplication
    subsection("Scalar multiplication (a * 3)");
    pp::vector s1 = a * 3;
    std::cout << "  a*3 = " << s1.to_string() << "\n";
    check(s1[0] == 3 && s1[1] == 6 && s1[2] == 9, "a*3 = {3,6,9}");

    subsection("Scalar multiplication (2 * a)");
    pp::vector s2 = 2 * a;
    std::cout << "  2*a = " << s2.to_string() << "\n";
    check(s2[0] == 2 && s2[1] == 4 && s2[2] == 6, "2*a = {2,4,6}");

    // Scalar division
    subsection("Scalar division (b / 2)");
    pp::vector s3 = b / 2;
    std::cout << "  b/2 = " << s3.to_string() << "\n";
    check(s3[0] == 2 && s3[1] == 2.5 && s3[2] == 3, "b/2 = {2,2.5,3}");

    // Compound assignment
    subsection("Compound assignment (+=, -=, *=, /=)");
    pp::vector ca = {1, 2, 3};
    ca += b;  check(ca[0] == 5 && ca[1] == 7 && ca[2] == 9, "a += b → {5,7,9}");
    ca -= b;  check(ca[0] == 1 && ca[1] == 2 && ca[2] == 3, "a -= b → {1,2,3}");
    ca *= 10; check(ca[0] == 10 && ca[1] == 20 && ca[2] == 30, "a *= 10 → {10,20,30}");
    ca /= 5;  check(ca[0] == 2 && ca[1] == 4 && ca[2] == 6, "a /= 5 → {2,4,6}");
}

// ═════════════════════════════════════════════════════════════════════
//  3. pp::vector — dot, norm, map, approx, print
// ═════════════════════════════════════════════════════════════════════
void test_vector_functions() {
    section("pp::vector — dot, norm, map, approx, print");

    pp::vector a = {1, 2, 3};
    pp::vector b = {4, 5, 6};

    // Dot product
    subsection("Dot product");
    NUMBER d1 = a.dot(b);         // member
    NUMBER d2 = pp::dot(a, b);    // free function
    std::cout << "    a·b (member) = " << d1 << "\n";
    std::cout << "    a·b (free)   = " << d2 << "\n";
    check(d1 == 32, "a.dot(b) = 1*4 + 2*5 + 3*6 = 32");
    check(d2 == 32, "dot(a,b) = 32 (free function)");

    // Norm
    subsection("Norm");
    pp::vector n = {3, 4};
    std::cout << "    ||{3,4}|| = " << n.norm() << "\n";
    check(std::abs(n.norm() - 5.0) < 1e-12, "||{3,4}|| = 5");

    // Map
    subsection("Map (apply function element-wise)");
    pp::vector sq = {1, 4, 9, 16};
    pp::vector sqr = sq.map([](NUMBER x){ return std::sqrt(x); });
    std::cout << "    sqrt({1,4,9,16}) = " << sqr.to_string() << "\n";
    check(std::abs(sqr[0]-1) < 1e-12 && std::abs(sqr[1]-2) < 1e-12 &&
          std::abs(sqr[2]-3) < 1e-12 && std::abs(sqr[3]-4) < 1e-12,
          "map(sqrt) = {1,2,3,4}");

    pp::vector dbl = a.map([](NUMBER x){ return 2*x; });
    std::cout << "    2*{1,2,3}        = " << dbl.to_string() << "\n";
    check(dbl[0] == 2 && dbl[1] == 4 && dbl[2] == 6, "map(2x) = {2,4,6}");

    // Approx
    subsection("Approximate equality (approx)");
    pp::vector c = {1.0 + 1e-8, 2.0 - 1e-8, 3.0 + 1e-8};
    check(pp::approx(a, a), "approx(a, a) = true (identical)");
    check(pp::approx(a, c), "approx(a, a+tiny) = true (within default tol)");
    check(!pp::approx(a, b), "approx(a, b) = false (different vectors)");
    pp::vector e = {1.1, 2.1, 3.1};
    check(!pp::approx(a, e), "approx off by 0.1 with default tol → false");
    check(pp::approx(a, e, 0.2, 0.0), "approx off by 0.1 with acc=0.2 → true");
    pp::vector f = {1, 2};
    check(!pp::approx(a, f), "approx with different sizes → false");

    // Print
    subsection("Print (visual check)");
    std::cout << "    v.print(\"vec = \"): ";
    a.print("vec = ");

    // to_string
    subsection("to_string");
    std::cout << "    a.to_string() = " << a.to_string() << "\n";
    check(a.to_string().size() > 0, "to_string produces non-empty output");
}

// ═════════════════════════════════════════════════════════════════════
//  4. pp::matrix — construction, element access, resize, setid
// ═════════════════════════════════════════════════════════════════════
void test_matrix_basics() {
    section("pp::matrix — basics");

    // Construction
    subsection("Construction");
    pp::matrix M0;
    check(M0.rows() == 0 && M0.columns() == 0, "Default constructor → 0×0");

    pp::matrix M(3, 4);
    check(M.rows() == 3 && M.columns() == 4, "matrix(3,4) → 3 rows, 4 cols");
    check(M[0, 0] == 0 && M[2, 3] == 0, "matrix(3,4) zero-initialized");

    // Element access: operator[i,j], operator(i,j), get/set
    subsection("Element access: operator[i,j], operator(i,j), get/set");
    M[1, 2] = 42;
    check(M[1, 2] == 42, "M[1,2] = 42 via operator[i,j]");
    check(M(1, 2) == 42, "M(1,2) = 42 via operator(i,j)");
    check(M.get(1, 2) == 42, "M.get(1,2) = 42");
    M.set(0, 3, -7.5);
    check(M.get(0, 3) == -7.5, "M.set(0,3,-7.5) then get(0,3) == -7.5");
    M(2, 0) = 99;
    check(M[2, 0] == 99, "M(2,0) = 99, read via M[2,0]");

    // Column access
    subsection("Column access (operator[j] returns column vector)");
    pp::matrix C(3, 2);
    C[0, 0] = 1; C[1, 0] = 2; C[2, 0] = 3;
    C[0, 1] = 4; C[1, 1] = 5; C[2, 1] = 6;
    pp::vector col0 = C[0]; // column 0
    check(col0[0] == 1 && col0[1] == 2 && col0[2] == 3, "Column 0 = {1,2,3}");

    // Resize
    subsection("Resize");
    pp::matrix R(2, 3);
    R[0,0]=1; R[0,1]=2; R[0,2]=3; R[1,0]=4; R[1,1]=5; R[1,2]=6;
    R.resize(4, 5);
    check(R.rows() == 4 && R.columns() == 5, "resize(4,5) → 4×5");
    check(R[0,0] == 1 && R[0,1] == 2 && R[1,2] == 6, "resize preserves existing data");
    check(R[2,0] == 0 && R[3,4] == 0, "resize zero-fills new elements");

    // setid
    subsection("Set to identity (setid)");
    pp::matrix I(4, 4);
    I[0,0] = 99; I[1,2] = 77;
    I.setid();
    check(is_identity(I), "setid turns arbitrary 4×4 into identity");
    I.print("    I = ");

    // Identity factory
    subsection("Identity factory (matrix::identity)");
    pp::matrix I3 = pp::matrix::identity(3);
    check(is_identity(I3), "matrix::identity(3) is 3×3 identity");

    // Random factory
    subsection("Random factory (matrix::random)");
    pp::matrix Rnd = pp::matrix::random(3, 4);
    check(Rnd.rows() == 3 && Rnd.columns() == 4, "random(3,4) has correct dimensions");
    bool has_nonzero = false;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            if (Rnd[i,j] != 0) has_nonzero = true;
    check(has_nonzero, "random matrix has non-zero entries");
}

// ═════════════════════════════════════════════════════════════════════
//  5. pp::matrix — arithmetic operations
// ═════════════════════════════════════════════════════════════════════
void test_matrix_arithmetic() {
    section("pp::matrix — arithmetic");

    pp::matrix A(2, 2); A[0,0]=1; A[0,1]=2; A[1,0]=3; A[1,1]=4;
    pp::matrix B(2, 2); B[0,0]=5; B[0,1]=6; B[1,0]=7; B[1,1]=8;

    std::cout << "    A = \n" << A.to_string();
    std::cout << "    B = \n" << B.to_string();

    // Addition
    subsection("Addition (A + B)");
    pp::matrix C = A + B;
    check(C[0,0]==6 && C[0,1]==8 && C[1,0]==10 && C[1,1]==12, "A+B correct");

    // Subtraction
    subsection("Subtraction (A - B)");
    pp::matrix D = A - B;
    check(D[0,0]==-4 && D[0,1]==-4 && D[1,0]==-4 && D[1,1]==-4, "A-B correct");

    // Unary negation
    subsection("Unary negation (-A)");
    pp::matrix N = -A;
    check(N[0,0]==-1 && N[0,1]==-2 && N[1,0]==-3 && N[1,1]==-4, "-A correct");

    // Scalar multiplication
    subsection("Scalar multiplication (A * 3,  2 * A)");
    pp::matrix S1 = A * 3;
    check(S1[0,0]==3 && S1[1,1]==12, "A*3 correct");
    pp::matrix S2 = 2 * A;
    check(S2[0,0]==2 && S2[1,1]==8, "2*A correct");

    // Scalar division
    subsection("Scalar division (B / 2)");
    pp::matrix S3 = B / 2;
    check(S3[0,0]==2.5 && S3[1,1]==4, "B/2 correct");

    // Compound assignment
    subsection("Compound assignment (+=, -=, *=, /=)");
    pp::matrix E = A;
    E += B; check(E[0,0]==6 && E[1,1]==12, "A += B");
    E -= B; check(E[0,0]==1 && E[1,1]==4, "A -= B back to original");
    E *= 10; check(E[0,0]==10 && E[1,1]==40, "A *= 10");
    E /= 5;  check(E[0,0]==2 && E[1,1]==8, "A /= 5");

    // Matrix-vector multiplication
    subsection("Matrix-vector multiplication (A * v)");
    pp::vector v = {1, 2};
    pp::vector r = A * v;
    std::cout << "    A * {1,2} = " << r.to_string() << "\n";
    check(r[0] == 5 && r[1] == 11,
          "A*v = {1*1+2*2, 3*1+4*2} = {5,11}");

    // Matrix-matrix multiplication
    subsection("Matrix-matrix multiplication (A * B)");
    pp::matrix AB = A * B;
    // [1 2][5 6] = [1*5+2*7  1*6+2*8] = [19 22]
    // [3 4][7 8]   [3*5+4*7  3*6+4*8]   [43 50]
    std::cout << "    A * B = \n" << AB.to_string();
    check(AB[0,0]==19 && AB[0,1]==22 && AB[1,0]==43 && AB[1,1]==50, "A*B correct");
}

// ═════════════════════════════════════════════════════════════════════
//  6. pp::matrix — transpose, T(), norm, print
// ═════════════════════════════════════════════════════════════════════
void test_matrix_functions() {
    section("pp::matrix — transpose, norm, print");

    pp::matrix A(2, 3);
    A[0,0]=1; A[0,1]=2; A[0,2]=3;
    A[1,0]=4; A[1,1]=5; A[1,2]=6;

    // Transpose
    subsection("Transpose");
    pp::matrix At = A.transpose();
    check(At.rows() == 3 && At.columns() == 2, "transpose swaps dimensions (2×3 → 3×2)");
    check(At[0,0]==1 && At[1,0]==2 && At[2,0]==3, "transpose row 0 → col 0");
    check(At[0,1]==4 && At[1,1]==5 && At[2,1]==6, "transpose row 1 → col 1");

    // T() alias
    subsection("T() (alias for transpose)");
    pp::matrix At2 = A.T();
    check(mat_approx(At, At2), "T() gives same result as transpose()");

    // Norm
    subsection("Frobenius norm");
    pp::matrix N(2, 2);
    N[0,0]=1; N[0,1]=2; N[1,0]=3; N[1,1]=4;
    // ||N||_F = sqrt(1+4+9+16) = sqrt(30)
    std::cout << "    ||[[1,2],[3,4]]||_F = " << N.norm() << "\n";
    check(std::abs(N.norm() - std::sqrt(30.0)) < 1e-12, "Frobenius norm = sqrt(30)");

    // Print
    subsection("Print (visual check)");
    A.print("    A = ");

    // to_string
    subsection("to_string");
    check(A.to_string().size() > 0, "to_string produces non-empty output");
}

// ═════════════════════════════════════════════════════════════════════
//  7. QR decomposition (modified Gram-Schmidt)
// ═════════════════════════════════════════════════════════════════════
void test_qr() {
    section("QR decomposition (modified Gram-Schmidt)");

    int n = 5, m = 5;
    pp::matrix A = pp::matrix::random(n, m);
    pp::qr qr(A);

    std::cout << "    A (" << n << "×" << m << ") random matrix\n\n";

    subsection("Decomposition properties");
    check(is_upper_tri(qr.R), "R is upper triangular");
    pp::matrix QtQ = qr.Q.transpose() * qr.Q;
    check(is_identity(QtQ), "Q^T Q = I (orthonormal columns)");
    check(mat_approx(qr.Q * qr.R, A), "Q R = A (reconstruction)");

    subsection("Linear solve (Ax = b)");
    pp::vector b = {1, 2, 3, 4, 5};
    pp::vector x = qr.solve(b);
    check(vec_approx(A * x, b), "A x = b");

    subsection("Determinant");
    NUMBER d = qr.det();
    std::cout << "    det(A) = " << d << "\n";
    check(std::abs(d) > 1e-15, "det is non-zero for random matrix");

    subsection("Inverse");
    pp::matrix Ainv = qr.inverse();
    check(is_identity(A * Ainv), "A * A^{-1} = I");

    // Tall matrix
    subsection("Tall matrix QR (8×4)");
    pp::matrix T = pp::matrix::random(8, 4);
    pp::qr qr_tall(T);
    check(is_upper_tri(qr_tall.R), "R is upper triangular (tall)");
    pp::matrix QtQ2 = qr_tall.Q.transpose() * qr_tall.Q;
    check(is_identity(QtQ2), "Q^T Q = I (tall)");
    check(mat_approx(qr_tall.Q * qr_tall.R, T), "Q R = A (tall)");
}

// ═════════════════════════════════════════════════════════════════════
//  8. LU decomposition (Doolittle with partial pivoting)
// ═════════════════════════════════════════════════════════════════════
void test_lu() {
    section("LU decomposition (Doolittle, partial pivoting)");

    int n = 5;
    pp::matrix A = pp::matrix::random(n, n);
    pp::lu lu(A);

    std::cout << "    A (" << n << "×" << n << ") random matrix\n\n";

    subsection("Decomposition properties");
    check(is_lower_tri(lu.L), "L is lower triangular");
    bool unit_diag = true;
    for (int i = 0; i < n; i++)
        if (std::abs(lu.L[i,i] - 1.0) > tol) { unit_diag = false; break; }
    check(unit_diag, "L has unit diagonal");
    check(is_upper_tri(lu.U), "U is upper triangular");

    // PA = LU
    pp::matrix PA(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            PA[i, j] = A[lu.perm[i], j];
    pp::matrix LU = lu.L * lu.U;
    check(mat_approx(PA, LU), "P A = L U");

    subsection("Linear solve (Ax = b)");
    pp::vector b = {1, 2, 3, 4, 5};
    pp::vector x = lu.solve(b);
    check(vec_approx(A * x, b), "A x = b");

    subsection("Determinant");
    NUMBER d = lu.det();
    std::cout << "    det(A) = " << d << "\n";
    check(std::abs(d) > 1e-15, "det is non-zero");

    subsection("Inverse");
    pp::matrix Ainv = lu.inverse();
    check(is_identity(A * Ainv), "A * A^{-1} = I");
}

// ═════════════════════════════════════════════════════════════════════
//  9. Cholesky decomposition (A = LL^T)
// ═════════════════════════════════════════════════════════════════════
void test_cholesky() {
    section("Cholesky decomposition (A = LL^T)");

    // Build SPD matrix: A = M^T M + I
    int n = 5;
    pp::matrix M = pp::matrix::random(n, n);
    pp::matrix A = M.transpose() * M + pp::matrix::identity(n);
    pp::cholesky chol(A);

    std::cout << "    A (" << n << "×" << n << ") symmetric positive definite\n\n";

    subsection("Decomposition properties");
    check(is_lower_tri(chol.L), "L is lower triangular");
    check(mat_approx(chol.L * chol.L.transpose(), A), "L L^T = A");

    subsection("Linear solve (Ax = b)");
    pp::vector b = {1, 2, 3, 4, 5};
    pp::vector x = chol.solve(b);
    check(vec_approx(A * x, b), "A x = b");

    subsection("Determinant");
    NUMBER d = chol.det();
    std::cout << "    det(A) = " << d << "\n";
    check(d > 0, "det > 0 for SPD matrix");

    subsection("Inverse");
    pp::matrix Ainv = chol.inverse();
    check(is_identity(A * Ainv), "A * A^{-1} = I");
}

// ═════════════════════════════════════════════════════════════════════
// 10. Householder QR decomposition
// ═════════════════════════════════════════════════════════════════════
void test_householder() {
    section("QR decomposition (Householder reflections)");

    int n = 5, m = 5;
    pp::matrix A = pp::matrix::random(n, m);
    pp::householder_qr hh(A);

    std::cout << "    A (" << n << "×" << m << ") random matrix\n\n";

    subsection("Decomposition properties");
    check(is_upper_tri(hh.R), "R is upper triangular");
    pp::matrix QtQ = hh.Q.transpose() * hh.Q;
    check(is_identity(QtQ), "Q^T Q = I");
    check(mat_approx(hh.Q * hh.R, A), "Q R = A");

    subsection("Linear solve (Ax = b)");
    pp::vector b = {1, 2, 3, 4, 5};
    pp::vector x = hh.solve(b);
    check(vec_approx(A * x, b), "A x = b");

    subsection("Determinant");
    NUMBER d = hh.det();
    std::cout << "    det(A) = " << d << "\n";
    check(std::abs(d) > 1e-15, "det is non-zero");

    subsection("Inverse");
    pp::matrix Ainv = hh.inverse();
    check(is_identity(A * Ainv), "A * A^{-1} = I");

    // Tall matrix
    subsection("Tall matrix QR (8×4)");
    pp::matrix T = pp::matrix::random(8, 4);
    pp::householder_qr hh2(T);
    check(is_upper_tri(hh2.R), "R is upper triangular (tall)");
    check(is_identity(hh2.Q.transpose() * hh2.Q), "Q^T Q = I (tall)");
    check(mat_approx(hh2.Q * hh2.R, T), "Q R = A (tall)");
}

// ═════════════════════════════════════════════════════════════════════
// 11. Givens QR decomposition
// ═════════════════════════════════════════════════════════════════════
void test_givens() {
    section("QR decomposition (Givens rotations)");

    int n = 5, m = 5;
    pp::matrix A = pp::matrix::random(n, m);
    pp::givens_qr gv(A);

    std::cout << "    A (" << n << "×" << m << ") random matrix\n\n";

    subsection("Decomposition properties");
    check(is_upper_tri(gv.R), "R is upper triangular");
    pp::matrix QtQ = gv.Q.transpose() * gv.Q;
    check(is_identity(QtQ), "Q^T Q = I");
    check(mat_approx(gv.Q * gv.R, A), "Q R = A");

    subsection("Linear solve (Ax = b)");
    pp::vector b = {1, 2, 3, 4, 5};
    pp::vector x = gv.solve(b);
    check(vec_approx(A * x, b), "A x = b");

    subsection("Determinant");
    NUMBER d = gv.det();
    std::cout << "    det(A) = " << d << "\n";
    check(std::abs(d) > 1e-15, "det is non-zero");

    subsection("Inverse");
    pp::matrix Ainv = gv.inverse();
    check(is_identity(A * Ainv), "A * A^{-1} = I");

    // Tall matrix
    subsection("Tall matrix QR (8×4)");
    pp::matrix T = pp::matrix::random(8, 4);
    pp::givens_qr gv2(T);
    check(is_upper_tri(gv2.R), "R is upper triangular (tall)");
    check(is_identity(gv2.Q.transpose() * gv2.Q), "Q^T Q = I (tall)");
    check(mat_approx(gv2.Q * gv2.R, T), "Q R = A (tall)");
}

// ═════════════════════════════════════════════════════════════════════
// 12. Cross-validation: determinant agreement across methods
// ═════════════════════════════════════════════════════════════════════
void test_determinant_cross() {
    section("Cross-validation: determinant agreement");

    pp::matrix A = pp::matrix::random(5, 5);

    pp::qr            qr(A);
    pp::lu             lu(A);
    pp::householder_qr hh(A);
    pp::givens_qr      gv(A);

    NUMBER d_qr = qr.det();
    NUMBER d_lu = lu.det();
    NUMBER d_hh = hh.det();
    NUMBER d_gv = gv.det();

    std::cout << "    det (QR Gram-Schmidt)       = " << d_qr << "\n";
    std::cout << "    det (LU partial pivoting)    = " << d_lu << "\n";
    std::cout << "    det (Householder)            = " << d_hh << "\n";
    std::cout << "    det (Givens)                 = " << d_gv << "\n\n";

    check(std::abs(d_qr - d_lu) < 1e-9, "QR det ≈ LU det");
    check(std::abs(d_qr - d_hh) < 1e-9, "QR det ≈ Householder det");
    check(std::abs(d_qr - d_gv) < 1e-9, "QR det ≈ Givens det");
}

// ═════════════════════════════════════════════════════════════════════
//  Main
// ═════════════════════════════════════════════════════════════════════
int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║     pp:: Linear Algebra Library — Full Test Suite       ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";

    // Vector tests
    test_vector_basics();
    test_vector_arithmetic();
    test_vector_functions();

    // Matrix tests
    test_matrix_basics();
    test_matrix_arithmetic();
    test_matrix_functions();

    // Decomposition tests
    test_qr();
    test_lu();
    test_cholesky();
    test_householder();
    test_givens();

    // Cross-validation
    test_determinant_cross();

    // Summary
    std::cout << "\n══════════════════════════════════════════════════════════\n";
    std::cout << "  TOTAL:  " << total_pass << " PASS,  " << total_fail << " FAIL\n";
    std::cout << "══════════════════════════════════════════════════════════\n";

    return total_fail;
}
