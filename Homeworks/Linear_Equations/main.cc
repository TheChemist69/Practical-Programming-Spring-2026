// Nikolaos Antonios Iakynthos Nemet
// This project was created using a combination of Claude Code and GitHub Copilot.

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <random>
#include "vector.h"
#include "matrix.h"
#include "qr.h"
#include "lu.h"
#include "cholesky.h"
#include "householder.h"
#include "givens.h"

static const double tol = 1e-9;

// Check whether a matrix is approximately upper triangular
bool is_upper_triangular(const pp::matrix& R) {
    for (int j = 0; j < R.columns(); j++) {
        for (int i = j + 1; i < R.rows(); i++) {
            if (std::abs(R[i, j]) > tol) return false;
        }
    }
    return true;
}

// Check whether a matrix is approximately the identity
bool is_identity(const pp::matrix& M) {
    for (int i = 0; i < M.rows(); i++) {
        for (int j = 0; j < M.columns(); j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            if (std::abs(M[i, j] - expected) > tol) return false;
        }
    }
    return true;
}

// Check whether two matrices are approximately equal
bool approx_equal(const pp::matrix& A, const pp::matrix& B) {
    return (A - B).norm() < tol;
}

// Check whether two vectors are approximately equal
bool approx_equal(const pp::vector& a, const pp::vector& b) {
    return (a - b).norm() < tol;
}

void task_1() {
    std::cout << "=== Task 1: QR-decomposition and linear solve ===\n\n";

    // Generate a random tall matrix A (n > m)
    int n = 8, m = 5;
    pp::matrix A = pp::matrix::random(n, m);

    std::cout << "Random tall matrix A (" << n << "x" << m << "):\n";
    std::cout << A.to_string() << "\n";

    // QR-decompose
    pp::qr decomp(A);

    // Check R is upper triangular
    std::cout << "R is upper triangular: "
              << (is_upper_triangular(decomp.R) ? "PASS" : "FAIL") << "\n";

    // Check Q^T Q = I
    pp::matrix QtQ = decomp.Q.transpose() * decomp.Q;
    std::cout << "Q^T Q = I:             "
              << (is_identity(QtQ) ? "PASS" : "FAIL") << "\n";

    // Check QR = A
    pp::matrix QR = decomp.Q * decomp.R;
    std::cout << "QR = A:                "
              << (approx_equal(QR, A) ? "PASS" : "FAIL") << "\n";

    std::cout << "\n--- Solving QRx = b ---\n\n";

    // Generate a random square matrix and vector for solving
    int sz = 6;
    pp::matrix A2 = pp::matrix::random(sz, sz);
    pp::vector b(sz);
    {
        static std::mt19937 gen(123);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < sz; i++) b[i] = dist(gen);
    }

    std::cout << "Random square matrix A (" << sz << "x" << sz << "):\n";
    std::cout << A2.to_string() << "\n";
    std::cout << "Random vector b: " << b.to_string() << "\n\n";

    pp::qr decomp2(A2);
    pp::vector x = decomp2.solve(b);

    std::cout << "Solution x:      " << x.to_string() << "\n";

    pp::vector Ax = A2 * x;
    std::cout << "A*x:             " << Ax.to_string() << "\n";
    std::cout << "b:               " << b.to_string() << "\n";
    std::cout << "Ax = b:          "
              << (approx_equal(Ax, b) ? "PASS" : "FAIL") << "\n";

    std::cout << "\n--- Determinant ---\n\n";
    std::cout << "det(A) = " << decomp2.det() << "\n";
}

void task_2() {
    std::cout << "\n=== Task 2: Matrix inverse ===\n\n";

    int sz = 6;
    pp::matrix A = pp::matrix::random(sz, sz);

    std::cout << "Random square matrix A (" << sz << "x" << sz << "):\n";
    std::cout << A.to_string() << "\n";

    pp::qr decomp(A);
    pp::matrix B = decomp.inverse();

    std::cout << "Inverse B = A^{-1}:\n";
    std::cout << B.to_string() << "\n";

    pp::matrix AB = A * B;
    pp::matrix I  = pp::matrix::identity(sz);

    std::cout << "A * B:\n" << AB.to_string() << "\n";
    std::cout << "AB = I: " << (approx_equal(AB, I) ? "PASS" : "FAIL") << "\n";
}

void task_3(int N) {
    // Just do a QR-decomposition of an NxN random matrix.
    // Timing is handled externally by the POSIX time utility.
    pp::matrix A = pp::matrix::random(N, N);
    pp::qr decomp(A);
    // Print something small so the computation is not optimized away
    std::cout << "N = " << N << ", det = " << decomp.det() << "\n";
}

// Check whether a matrix is approximately lower triangular
bool is_lower_triangular(const pp::matrix& L) {
    for (int j = 0; j < L.columns(); j++) {
        for (int i = 0; i < j; i++) {
            if (std::abs(L[i, j]) > tol) return false;
        }
    }
    return true;
}

void task_4() {
    std::cout << "\n=== Task 4: LU-decomposition (Doolittle with partial pivoting) ===\n\n";

    int sz = 6;
    pp::matrix A = pp::matrix::random(sz, sz);

    std::cout << "Random square matrix A (" << sz << "x" << sz << "):\n";
    std::cout << A.to_string() << "\n";

    pp::lu decomp(A);

    // Check L is lower triangular with unit diagonal
    std::cout << "L is lower triangular: "
              << (is_lower_triangular(decomp.L) ? "PASS" : "FAIL") << "\n";

    bool unit_diag = true;
    for (int i = 0; i < sz; i++) {
        if (std::abs(decomp.L[i, i] - 1.0) > tol) { unit_diag = false; break; }
    }
    std::cout << "L has unit diagonal:   " << (unit_diag ? "PASS" : "FAIL") << "\n";

    // Check U is upper triangular
    std::cout << "U is upper triangular: "
              << (is_upper_triangular(decomp.U) ? "PASS" : "FAIL") << "\n";

    // Check PA = LU
    pp::matrix LU = decomp.L * decomp.U;
    pp::matrix PA(sz, sz);
    for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz; j++) {
            PA[i, j] = A[decomp.perm[i], j];
        }
    }
    std::cout << "PA = LU:               "
              << (approx_equal(PA, LU) ? "PASS" : "FAIL") << "\n";

    // Solve Ax = b
    pp::vector b(sz);
    {
        static std::mt19937 gen(456);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < sz; i++) b[i] = dist(gen);
    }
    pp::vector x = decomp.solve(b);
    pp::vector Ax = A * x;
    std::cout << "Ax = b:                "
              << (approx_equal(Ax, b) ? "PASS" : "FAIL") << "\n";

    // Determinant
    std::cout << "det(A) = " << decomp.det() << "\n";

    // Inverse
    pp::matrix B = decomp.inverse();
    pp::matrix AB = A * B;
    pp::matrix I = pp::matrix::identity(sz);
    std::cout << "AB = I:                "
              << (approx_equal(AB, I) ? "PASS" : "FAIL") << "\n";
}

void task_5() {
    std::cout << "\n=== Task 5: Cholesky decomposition ===\n\n";

    // Generate a symmetric positive definite matrix: A = M^T * M + I
    int sz = 6;
    pp::matrix M = pp::matrix::random(sz, sz);
    pp::matrix A = M.transpose() * M + pp::matrix::identity(sz);

    std::cout << "Symmetric positive definite matrix A (" << sz << "x" << sz << "):\n";
    std::cout << A.to_string() << "\n";

    pp::cholesky decomp(A);

    // Check L is lower triangular
    std::cout << "L is lower triangular: "
              << (is_lower_triangular(decomp.L) ? "PASS" : "FAIL") << "\n";

    // Check LL^T = A
    pp::matrix LLt = decomp.L * decomp.L.transpose();
    std::cout << "LL^T = A:              "
              << (approx_equal(LLt, A) ? "PASS" : "FAIL") << "\n";

    // Solve Ax = b
    pp::vector b(sz);
    {
        static std::mt19937 gen(789);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < sz; i++) b[i] = dist(gen);
    }
    pp::vector x = decomp.solve(b);
    pp::vector Ax_vec = A * x;
    std::cout << "Ax = b:                "
              << (approx_equal(Ax_vec, b) ? "PASS" : "FAIL") << "\n";

    // Determinant
    std::cout << "det(A) = " << decomp.det() << "\n";

    // Inverse
    pp::matrix B = decomp.inverse();
    pp::matrix AB = A * B;
    pp::matrix I = pp::matrix::identity(sz);
    std::cout << "AB = I:                "
              << (approx_equal(AB, I) ? "PASS" : "FAIL") << "\n";
}

void task_6() {
    std::cout << "\n=== Task 6: QR-decomposition via Householder reflections ===\n\n";

    // Tall matrix
    int n = 8, m = 5;
    pp::matrix A = pp::matrix::random(n, m);

    pp::householder_qr decomp(A);

    std::cout << "R is upper triangular: "
              << (is_upper_triangular(decomp.R) ? "PASS" : "FAIL") << "\n";

    pp::matrix QtQ = decomp.Q.transpose() * decomp.Q;
    std::cout << "Q^T Q = I:             "
              << (is_identity(QtQ) ? "PASS" : "FAIL") << "\n";

    pp::matrix QR = decomp.Q * decomp.R;
    std::cout << "QR = A:                "
              << (approx_equal(QR, A) ? "PASS" : "FAIL") << "\n";

    // Solve with a square matrix
    int sz = 6;
    pp::matrix A2 = pp::matrix::random(sz, sz);
    pp::vector b(sz);
    {
        static std::mt19937 gen(321);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < sz; i++) b[i] = dist(gen);
    }
    pp::householder_qr decomp2(A2);
    pp::vector x = decomp2.solve(b);
    pp::vector Ax = A2 * x;
    std::cout << "Ax = b:                "
              << (approx_equal(Ax, b) ? "PASS" : "FAIL") << "\n";

    // Inverse
    pp::matrix B = decomp2.inverse();
    pp::matrix AB = A2 * B;
    pp::matrix I = pp::matrix::identity(sz);
    std::cout << "AB = I:                "
              << (approx_equal(AB, I) ? "PASS" : "FAIL") << "\n";
}

void task_7() {
    std::cout << "\n=== Task 7: QR-decomposition via Givens rotations ===\n\n";

    // Tall matrix
    int n = 8, m = 5;
    pp::matrix A = pp::matrix::random(n, m);

    pp::givens_qr decomp(A);

    std::cout << "R is upper triangular: "
              << (is_upper_triangular(decomp.R) ? "PASS" : "FAIL") << "\n";

    pp::matrix QtQ = decomp.Q.transpose() * decomp.Q;
    std::cout << "Q^T Q = I:             "
              << (is_identity(QtQ) ? "PASS" : "FAIL") << "\n";

    pp::matrix QR = decomp.Q * decomp.R;
    std::cout << "QR = A:                "
              << (approx_equal(QR, A) ? "PASS" : "FAIL") << "\n";

    // Solve with a square matrix
    int sz = 6;
    pp::matrix A2 = pp::matrix::random(sz, sz);
    pp::vector b(sz);
    {
        static std::mt19937 gen(654);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < sz; i++) b[i] = dist(gen);
    }
    pp::givens_qr decomp2(A2);
    pp::vector x = decomp2.solve(b);
    pp::vector Ax = A2 * x;
    std::cout << "Ax = b:                "
              << (approx_equal(Ax, b) ? "PASS" : "FAIL") << "\n";

    // Inverse
    pp::matrix B = decomp2.inverse();
    pp::matrix AB = A2 * B;
    pp::matrix I = pp::matrix::identity(sz);
    std::cout << "AB = I:                "
              << (approx_equal(AB, I) ? "PASS" : "FAIL") << "\n";
}

int main(int argc, char* argv[]) {
    int size = 0; // 0 means run the demo/checks

    for (int i = 0; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-size" && i + 1 < argc) {
            size = std::stoi(argv[i + 1]);
        }
    }

    if (size > 0) {
        // Timing mode: just factorize an NxN matrix
        task_3(size);
    } else {
        // Demo mode: run all checks
        std::cout << std::fixed << std::setprecision(6);
        task_1();
        task_2();
        task_4();
        task_5();
        task_6();
        task_7();
    }

    return 0;
}
