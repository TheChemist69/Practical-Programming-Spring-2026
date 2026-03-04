// Comprehensive stress test for all decompositions.
// Tests edge cases, known matrices, determinant signs, and varied sizes.
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <string>
#include "vector.h"
#include "matrix.h"
#include "qr.h"
#include "lu.h"
#include "cholesky.h"
#include "householder.h"
#include "givens.h"

static const double tol = 1e-9;
static int pass_count = 0;
static int fail_count = 0;

void check(const std::string& label, bool ok) {
    std::cout << label << ": " << (ok ? "PASS" : "** FAIL **") << "\n";
    if (ok) pass_count++; else fail_count++;
}

bool approx_eq(double a, double b, double eps = tol) {
    return std::abs(a - b) < eps;
}

bool mat_approx_eq(const pp::matrix& A, const pp::matrix& B, double eps = tol) {
    if (A.rows() != B.rows() || A.columns() != B.columns()) return false;
    return (A - B).norm() < eps;
}

bool vec_approx_eq(const pp::vector& a, const pp::vector& b, double eps = tol) {
    return (a - b).norm() < eps;
}

bool is_identity(const pp::matrix& M, double eps = tol) {
    for (int i = 0; i < M.rows(); i++)
        for (int j = 0; j < M.columns(); j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            if (std::abs(M[i, j] - expected) > eps) return false;
        }
    return true;
}

// ─────── Test 1: 1x1 matrix ───────
void test_1x1() {
    std::cout << "\n=== Test: 1x1 matrix ===\n";
    pp::matrix A(1, 1);
    A[0, 0] = 5.0;
    pp::vector b({3.0});

    // QR (Gram-Schmidt)
    pp::qr qr_d(A);
    check("  QR solve 1x1", vec_approx_eq(A * qr_d.solve(b), b));
    check("  QR det 1x1", approx_eq(qr_d.det(), 5.0));
    check("  QR inv 1x1", approx_eq(qr_d.inverse()[0, 0], 0.2));

    // LU
    pp::lu lu_d(A);
    check("  LU solve 1x1", vec_approx_eq(A * lu_d.solve(b), b));
    check("  LU det 1x1", approx_eq(lu_d.det(), 5.0));
    check("  LU inv 1x1", approx_eq(lu_d.inverse()[0, 0], 0.2));

    // Cholesky (5 > 0, so SPD)
    pp::cholesky ch_d(A);
    check("  Chol solve 1x1", vec_approx_eq(A * ch_d.solve(b), b));
    check("  Chol det 1x1", approx_eq(ch_d.det(), 5.0));

    // Householder
    pp::householder_qr hh_d(A);
    check("  HH solve 1x1", vec_approx_eq(A * hh_d.solve(b), b));

    // Givens
    pp::givens_qr gv_d(A);
    check("  Giv solve 1x1", vec_approx_eq(A * gv_d.solve(b), b));
}

// ─────── Test 2: 2x2 known matrix ───────
void test_2x2_known() {
    std::cout << "\n=== Test: 2x2 known matrix [[1,2],[3,4]] ===\n";
    pp::matrix A(2, 2);
    A[0,0] = 1; A[0,1] = 2;
    A[1,0] = 3; A[1,1] = 4;
    // det = 1*4 - 2*3 = -2
    // inv = [[-2, 1], [1.5, -0.5]]
    pp::vector b({5.0, 11.0}); // x = [1, 2]

    // QR (Gram-Schmidt)
    pp::qr qr_d(A);
    pp::vector x_qr = qr_d.solve(b);
    check("  QR solve [[1,2],[3,4]]", vec_approx_eq(x_qr, pp::vector({1.0, 2.0})));
    // GS det is always positive (norms), so it returns |det| = 2
    double qr_det = qr_d.det();
    check("  QR |det| = 2", approx_eq(std::abs(qr_det), 2.0));
    std::cout << "    QR  det() returned: " << qr_det
              << " (expected -2, got " << qr_det << ")\n";

    // LU
    pp::lu lu_d(A);
    pp::vector x_lu = lu_d.solve(b);
    check("  LU solve [[1,2],[3,4]]", vec_approx_eq(x_lu, pp::vector({1.0, 2.0})));
    double lu_det = lu_d.det();
    check("  LU det = -2", approx_eq(lu_det, -2.0));
    std::cout << "    LU  det() returned: " << lu_det << "\n";

    // Inverse check
    pp::matrix inv_expected(2, 2);
    inv_expected[0,0] = -2.0; inv_expected[0,1] = 1.0;
    inv_expected[1,0] = 1.5;  inv_expected[1,1] = -0.5;
    pp::matrix inv_lu = lu_d.inverse();
    check("  LU inv [[1,2],[3,4]]", mat_approx_eq(inv_lu, inv_expected));

    // Householder
    pp::householder_qr hh_d(A);
    pp::vector x_hh = hh_d.solve(b);
    check("  HH solve [[1,2],[3,4]]", vec_approx_eq(x_hh, pp::vector({1.0, 2.0})));
    double hh_det = hh_d.det();
    check("  HH |det| = 2", approx_eq(std::abs(hh_det), 2.0));
    std::cout << "    HH  det() returned: " << hh_det
              << " (expected -2, got " << hh_det << ")\n";

    // Givens
    pp::givens_qr gv_d(A);
    pp::vector x_gv = gv_d.solve(b);
    check("  Giv solve [[1,2],[3,4]]", vec_approx_eq(x_gv, pp::vector({1.0, 2.0})));
    double gv_det = gv_d.det();
    check("  Giv |det| = 2", approx_eq(std::abs(gv_det), 2.0));
    std::cout << "    Giv det() returned: " << gv_det
              << " (expected -2, got " << gv_det << ")\n";
}

// ─────── Test 3: Identity matrix ───────
void test_identity() {
    std::cout << "\n=== Test: Identity matrix (N=5) ===\n";
    int n = 5;
    pp::matrix I = pp::matrix::identity(n);
    pp::vector b({1, 2, 3, 4, 5});

    pp::qr qr_d(I);
    check("  QR solve I*x=b", vec_approx_eq(qr_d.solve(b), b));
    check("  QR det(I) = 1", approx_eq(qr_d.det(), 1.0));
    check("  QR inv(I) = I", is_identity(qr_d.inverse()));

    pp::lu lu_d(I);
    check("  LU solve I*x=b", vec_approx_eq(lu_d.solve(b), b));
    check("  LU det(I) = 1", approx_eq(lu_d.det(), 1.0));
    check("  LU inv(I) = I", is_identity(lu_d.inverse()));

    pp::cholesky ch_d(I);
    check("  Chol solve I*x=b", vec_approx_eq(ch_d.solve(b), b));
    check("  Chol det(I) = 1", approx_eq(ch_d.det(), 1.0));
    check("  Chol inv(I) = I", is_identity(ch_d.inverse()));

    pp::householder_qr hh_d(I);
    check("  HH solve I*x=b", vec_approx_eq(hh_d.solve(b), b));
    check("  HH inv(I) = I", is_identity(hh_d.inverse()));

    pp::givens_qr gv_d(I);
    check("  Giv solve I*x=b", vec_approx_eq(gv_d.solve(b), b));
    check("  Giv inv(I) = I", is_identity(gv_d.inverse()));
}

// ─────── Test 4: Diagonal matrix ───────
void test_diagonal() {
    std::cout << "\n=== Test: Diagonal matrix diag(2,3,4) ===\n";
    int n = 3;
    pp::matrix D(n, n);
    D[0,0] = 2; D[1,1] = 3; D[2,2] = 4;
    pp::vector b({6, 9, 12}); // x = {3, 3, 3}
    pp::vector x_exp({3, 3, 3});

    pp::qr qr_d(D);
    check("  QR solve diag", vec_approx_eq(qr_d.solve(b), x_exp));
    check("  QR det diag = 24", approx_eq(qr_d.det(), 24.0));

    pp::lu lu_d(D);
    check("  LU solve diag", vec_approx_eq(lu_d.solve(b), x_exp));
    check("  LU det diag = 24", approx_eq(lu_d.det(), 24.0));

    pp::cholesky ch_d(D);
    check("  Chol solve diag", vec_approx_eq(ch_d.solve(b), x_exp));
    check("  Chol det diag = 24", approx_eq(ch_d.det(), 24.0));

    pp::householder_qr hh_d(D);
    check("  HH solve diag", vec_approx_eq(hh_d.solve(b), x_exp));

    pp::givens_qr gv_d(D);
    check("  Giv solve diag", vec_approx_eq(gv_d.solve(b), x_exp));
}

// ─────── Test 5: Permutation matrix (zeros on diagonal) ───────
void test_permutation_matrix() {
    std::cout << "\n=== Test: Permutation matrix (all-zero diagonal) ===\n";
    // P = [[0,1,0],[0,0,1],[1,0,0]], det = 1
    pp::matrix P(3, 3);
    P[0,1] = 1; P[1,2] = 1; P[2,0] = 1;
    pp::vector b({7, 8, 9}); // Px = b => x = P^T b = {9, 7, 8}
    pp::vector x_exp({9, 7, 8});

    pp::lu lu_d(P);
    check("  LU solve permutation", vec_approx_eq(lu_d.solve(b), x_exp));
    check("  LU det permutation = 1", approx_eq(lu_d.det(), 1.0));

    pp::qr qr_d(P);
    check("  QR solve permutation", vec_approx_eq(qr_d.solve(b), x_exp));

    pp::householder_qr hh_d(P);
    check("  HH solve permutation", vec_approx_eq(hh_d.solve(b), x_exp));

    pp::givens_qr gv_d(P);
    check("  Giv solve permutation", vec_approx_eq(gv_d.solve(b), x_exp));
}

// ─────── Test 6: Negative determinant ───────
void test_negative_det() {
    std::cout << "\n=== Test: Determinant sign ===\n";
    // A = [[0,1],[1,0]], det = -1
    pp::matrix A(2, 2);
    A[0,1] = 1; A[1,0] = 1;

    pp::lu lu_d(A);
    double lu_det = lu_d.det();
    check("  LU det [[0,1],[1,0]] = -1", approx_eq(lu_det, -1.0));
    std::cout << "    LU  det = " << lu_det << "\n";

    pp::qr qr_d(A);
    double qr_det = qr_d.det();
    std::cout << "    QR  det = " << qr_det << " (should be -1, GS always positive?)\n";
    // GS det returns product of R diag which are norms → always positive
    // So this should return 1, not -1 — a known limitation
    bool qr_sign_correct = approx_eq(qr_det, -1.0);
    bool qr_abs_correct = approx_eq(std::abs(qr_det), 1.0);
    if (!qr_sign_correct) {
        std::cout << "    ** QR det has WRONG SIGN (returns |det| not det) **\n";
        fail_count++;
    } else {
        pass_count++;
    }
    check("  QR |det| correct", qr_abs_correct);

    pp::householder_qr hh_d(A);
    double hh_det = hh_d.det();
    std::cout << "    HH  det = " << hh_det << " (should be -1)\n";
    bool hh_sign_correct = approx_eq(hh_det, -1.0);
    if (!hh_sign_correct) {
        std::cout << "    ** HH det has WRONG SIGN **\n";
        fail_count++;
    } else {
        pass_count++;
    }

    pp::givens_qr gv_d(A);
    double gv_det = gv_d.det();
    std::cout << "    Giv det = " << gv_det << " (should be -1)\n";
    bool gv_sign_correct = approx_eq(gv_det, -1.0);
    if (!gv_sign_correct) {
        std::cout << "    ** Giv det has WRONG SIGN **\n";
        fail_count++;
    } else {
        pass_count++;
    }
}

// ─────── Test 7: Larger random matrices (N=20, 50) ───────
void test_larger_random() {
    std::cout << "\n=== Test: Larger random matrices ===\n";
    std::mt19937 gen(999);
    std::uniform_real_distribution<double> dist(-5.0, 5.0);

    for (int n : {10, 20, 50}) {
        // Generate random matrix
        pp::matrix A(n, n);
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
                A[i, j] = dist(gen);

        // Generate random b
        pp::vector b(n);
        for (int i = 0; i < n; i++) b[i] = dist(gen);

        std::string prefix = "  N=" + std::to_string(n);

        // QR
        pp::qr qr_d(A);
        check(prefix + " QR Ax=b", vec_approx_eq(A * qr_d.solve(b), b));
        check(prefix + " QR A*inv=I", is_identity(A * qr_d.inverse()));

        // LU
        pp::lu lu_d(A);
        check(prefix + " LU Ax=b", vec_approx_eq(A * lu_d.solve(b), b));
        check(prefix + " LU A*inv=I", is_identity(A * lu_d.inverse()));

        // Householder
        pp::householder_qr hh_d(A);
        check(prefix + " HH Ax=b", vec_approx_eq(A * hh_d.solve(b), b));
        check(prefix + " HH A*inv=I", is_identity(A * hh_d.inverse()));

        // Givens
        pp::givens_qr gv_d(A);
        check(prefix + " Giv Ax=b", vec_approx_eq(A * gv_d.solve(b), b));
        check(prefix + " Giv A*inv=I", is_identity(A * gv_d.inverse()));

        // Cholesky (need SPD: A^T A + I)
        pp::matrix S = A.transpose() * A + pp::matrix::identity(n);
        pp::cholesky ch_d(S);
        check(prefix + " Chol Sx=b", vec_approx_eq(S * ch_d.solve(b), b));
        check(prefix + " Chol S*inv=I", is_identity(S * ch_d.inverse()));

        // Cross-check: LU det vs QR |det|
        double lu_det_val = lu_d.det();
        double qr_det_val = qr_d.det();
        check(prefix + " |LU det| == |QR det|",
              approx_eq(std::abs(lu_det_val), std::abs(qr_det_val),
                        std::max(1e-6 * std::abs(lu_det_val), 1e-6)));
    }
}

// ─────── Test 8: Hilbert matrix (ill-conditioned) ───────
void test_hilbert() {
    std::cout << "\n=== Test: Hilbert matrix (ill-conditioned, N=5) ===\n";
    int n = 5;
    pp::matrix H(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            H[i, j] = 1.0 / (i + j + 1);

    // b = H * [1,1,...,1] so x = [1,1,...,1]
    pp::vector ones(n);
    for (int i = 0; i < n; i++) ones[i] = 1.0;
    pp::vector b = H * ones;

    double loose_tol = 1e-4; // Hilbert is ill-conditioned

    pp::lu lu_d(H);
    check("  LU solve Hilbert", vec_approx_eq(lu_d.solve(b), ones, loose_tol));

    pp::qr qr_d(H);
    check("  QR solve Hilbert", vec_approx_eq(qr_d.solve(b), ones, loose_tol));

    pp::householder_qr hh_d(H);
    check("  HH solve Hilbert", vec_approx_eq(hh_d.solve(b), ones, loose_tol));

    pp::givens_qr gv_d(H);
    check("  Giv solve Hilbert", vec_approx_eq(gv_d.solve(b), ones, loose_tol));

    // Hilbert is SPD
    pp::cholesky ch_d(H);
    check("  Chol solve Hilbert", vec_approx_eq(ch_d.solve(b), ones, loose_tol));
}

// ─────── Test 9: Tall matrix QR (overdetermined least squares) ───────
void test_tall_qr() {
    std::cout << "\n=== Test: Tall matrix QR (10x3) ===\n";
    int n = 10, m = 3;
    std::mt19937 gen(777);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    pp::matrix A(n, m);
    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++)
            A[i, j] = dist(gen);

    pp::qr qr_d(A);
    pp::matrix QR = qr_d.Q * qr_d.R;
    check("  QR=A for tall", mat_approx_eq(QR, A));
    pp::matrix QtQ = qr_d.Q.transpose() * qr_d.Q;
    check("  Q^TQ=I for tall", is_identity(QtQ));

    pp::householder_qr hh_d(A);
    pp::matrix QR_hh = hh_d.Q * hh_d.R;
    check("  HH QR=A for tall", mat_approx_eq(QR_hh, A));
    pp::matrix QtQ_hh = hh_d.Q.transpose() * hh_d.Q;
    check("  HH Q^TQ=I for tall", is_identity(QtQ_hh));

    pp::givens_qr gv_d(A);
    pp::matrix QR_gv = gv_d.Q * gv_d.R;
    check("  Giv QR=A for tall", mat_approx_eq(QR_gv, A));
    pp::matrix QtQ_gv = gv_d.Q.transpose() * gv_d.Q;
    check("  Giv Q^TQ=I for tall", is_identity(QtQ_gv));
}

// ─────── Test 10: Vector operations ───────
void test_vector_ops() {
    std::cout << "\n=== Test: Vector operations ===\n";
    pp::vector a({1, 2, 3});
    pp::vector b({4, 5, 6});

    check("  a+b", vec_approx_eq(a + b, pp::vector({5, 7, 9})));
    check("  a-b", vec_approx_eq(a - b, pp::vector({-3, -3, -3})));
    check("  a*2", vec_approx_eq(a * 2, pp::vector({2, 4, 6})));
    check("  2*a", vec_approx_eq(2.0 * a, pp::vector({2, 4, 6})));
    check("  a/2", vec_approx_eq(a / 2.0, pp::vector({0.5, 1.0, 1.5})));
    check("  -a", vec_approx_eq(-a, pp::vector({-1, -2, -3})));
    check("  dot(a,b)=32", approx_eq(a.dot(b), 32.0));
    check("  norm({3,4})=5", approx_eq(pp::vector({3.0, 4.0}).norm(), 5.0));

    // Compound assignment
    pp::vector c = a;
    c += b;
    check("  a+=b", vec_approx_eq(c, pp::vector({5, 7, 9})));
    c = a;
    c -= b;
    check("  a-=b", vec_approx_eq(c, pp::vector({-3, -3, -3})));
    c = a;
    c *= 3.0;
    check("  a*=3", vec_approx_eq(c, pp::vector({3, 6, 9})));
    c = a;
    c /= 2.0;
    check("  a/=2", vec_approx_eq(c, pp::vector({0.5, 1.0, 1.5})));
}

// ─────── Test 11: Matrix operations ───────
void test_matrix_ops() {
    std::cout << "\n=== Test: Matrix operations ===\n";
    pp::matrix A(2, 2);
    A[0,0]=1; A[0,1]=2; A[1,0]=3; A[1,1]=4;

    pp::matrix B(2, 2);
    B[0,0]=5; B[0,1]=6; B[1,0]=7; B[1,1]=8;

    // A+B
    pp::matrix C = A + B;
    check("  A+B [0,0]=6", approx_eq(C[0,0], 6.0));
    check("  A+B [1,1]=12", approx_eq(C[1,1], 12.0));

    // A*B = [[19,22],[43,50]]
    pp::matrix D = A * B;
    check("  A*B [0,0]=19", approx_eq(D[0,0], 19.0));
    check("  A*B [0,1]=22", approx_eq(D[0,1], 22.0));
    check("  A*B [1,0]=43", approx_eq(D[1,0], 43.0));
    check("  A*B [1,1]=50", approx_eq(D[1,1], 50.0));

    // Transpose
    pp::matrix At = A.transpose();
    check("  A^T [0,1]=3", approx_eq(At[0,1], 3.0));
    check("  A^T [1,0]=2", approx_eq(At[1,0], 2.0));

    // Non-square transpose
    pp::matrix R(2, 3);
    R[0,0]=1; R[0,1]=2; R[0,2]=3; R[1,0]=4; R[1,1]=5; R[1,2]=6;
    pp::matrix Rt = R.transpose();
    check("  R^T shape 3x2", Rt.rows() == 3 && Rt.columns() == 2);
    check("  R^T [2,0]=3", approx_eq(Rt[2,0], 3.0));
    check("  R^T [2,1]=6", approx_eq(Rt[2,1], 6.0));

    // Frobenius norm
    // A = [[1,2],[3,4]], ||A||_F = sqrt(1+4+9+16) = sqrt(30)
    check("  Frobenius norm", approx_eq(A.norm(), std::sqrt(30.0)));

    // Scalar ops
    pp::matrix E = A * 2.0;
    check("  A*2 [0,0]=2", approx_eq(E[0,0], 2.0));
    check("  A*2 [1,1]=8", approx_eq(E[1,1], 8.0));

    pp::matrix F = A / 2.0;
    check("  A/2 [0,0]=0.5", approx_eq(F[0,0], 0.5));
}

// ─────── Test 12: Cross-validate determinants across decompositions ───────
void test_det_cross_validation() {
    std::cout << "\n=== Test: Determinant cross-validation ===\n";
    std::mt19937 gen(12345);
    std::uniform_real_distribution<double> dist(-2.0, 2.0);

    for (int trial = 0; trial < 10; trial++) {
        int n = 4 + trial; // sizes 4..13
        pp::matrix A(n, n);
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
                A[i, j] = dist(gen);

        pp::lu lu_d(A);
        pp::qr qr_d(A);
        pp::householder_qr hh_d(A);
        pp::givens_qr gv_d(A);

        double lu_det = lu_d.det();
        double qr_det = qr_d.det();
        double hh_det = hh_d.det();
        double gv_det = gv_d.det();

        double eps = std::max(1e-6 * std::abs(lu_det), 1e-9);
        std::string pfx = "  N=" + std::to_string(n);

        // LU should give correct signed det
        // QR (GS) gives |det| (product of norms, always positive)
        // Householder gives signed det(R), sign may differ
        // Givens gives det(R), and det(Q_givens)=1, so should be correct

        check(pfx + " |QR|==|LU|", approx_eq(std::abs(qr_det), std::abs(lu_det), eps));
        check(pfx + " |HH|==|LU|", approx_eq(std::abs(hh_det), std::abs(lu_det), eps));
        check(pfx + " Giv==LU", approx_eq(gv_det, lu_det, eps));

        // Check if QR sign matches LU
        if (!approx_eq(qr_det, lu_det, eps)) {
            std::cout << "      QR sign differs: QR=" << qr_det
                      << " LU=" << lu_det << "\n";
        }
        // Check if HH sign matches LU
        if (!approx_eq(hh_det, lu_det, eps)) {
            std::cout << "      HH sign differs: HH=" << hh_det
                      << " LU=" << lu_det << "\n";
        }
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(9);

    test_vector_ops();
    test_matrix_ops();
    test_1x1();
    test_2x2_known();
    test_identity();
    test_diagonal();
    test_permutation_matrix();
    test_negative_det();
    test_larger_random();
    test_hilbert();
    test_tall_qr();
    test_det_cross_validation();

    std::cout << "\n════════════════════════════════════════\n";
    std::cout << "  TOTAL PASS: " << pass_count << "\n";
    std::cout << "  TOTAL FAIL: " << fail_count << "\n";
    std::cout << "════════════════════════════════════════\n";

    return fail_count > 0 ? 1 : 0;
}
