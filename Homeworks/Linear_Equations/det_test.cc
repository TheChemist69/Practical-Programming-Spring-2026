// Quick determinant sign verification
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include "qr.h"
#include "lu.h"
#include "householder.h"
#include "givens.h"

int main() {
    std::cout << std::fixed << std::setprecision(9);
    int fails = 0;

    // Test 1: [[0,1],[1,0]], det = -1
    {
        pp::matrix A(2,2);
        A[0,1]=1; A[1,0]=1;
        double lu_d = pp::lu(A).det();
        double qr_d = pp::qr(A).det();
        double hh_d = pp::householder_qr(A).det();
        double gv_d = pp::givens_qr(A).det();
        std::cout << "[[0,1],[1,0]] det=-1:  LU=" << lu_d
                  << " QR=" << qr_d << " HH=" << hh_d << " Giv=" << gv_d << "\n";
        if (std::abs(qr_d - (-1.0)) > 1e-9) { std::cout << "  QR FAIL\n"; fails++; }
        if (std::abs(hh_d - (-1.0)) > 1e-9) { std::cout << "  HH FAIL\n"; fails++; }
    }

    // Test 2: [[1,2],[3,4]], det = -2
    {
        pp::matrix A(2,2);
        A[0,0]=1; A[0,1]=2; A[1,0]=3; A[1,1]=4;
        double lu_d = pp::lu(A).det();
        double qr_d = pp::qr(A).det();
        double hh_d = pp::householder_qr(A).det();
        double gv_d = pp::givens_qr(A).det();
        std::cout << "[[1,2],[3,4]] det=-2:  LU=" << lu_d
                  << " QR=" << qr_d << " HH=" << hh_d << " Giv=" << gv_d << "\n";
        if (std::abs(qr_d - (-2.0)) > 1e-9) { std::cout << "  QR FAIL\n"; fails++; }
        if (std::abs(hh_d - (-2.0)) > 1e-9) { std::cout << "  HH FAIL\n"; fails++; }
    }

    // Test 3: Cross-validate on random matrices (N=4..13)
    std::mt19937 gen(12345);
    std::uniform_real_distribution<double> dist(-2.0, 2.0);
    for (int n = 4; n <= 13; n++) {
        pp::matrix A(n, n);
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
                A[i,j] = dist(gen);

        double lu_d = pp::lu(A).det();
        double qr_d = pp::qr(A).det();
        double hh_d = pp::householder_qr(A).det();
        double gv_d = pp::givens_qr(A).det();
        double eps = std::max(1e-6 * std::abs(lu_d), 1e-9);

        bool qr_ok = std::abs(qr_d - lu_d) < eps;
        bool hh_ok = std::abs(hh_d - lu_d) < eps;
        bool gv_ok = std::abs(gv_d - lu_d) < eps;

        std::cout << "N=" << n << ": LU=" << lu_d << " QR=" << qr_d
                  << " HH=" << hh_d << " Giv=" << gv_d;
        if (!qr_ok) { std::cout << " QR-MISMATCH"; fails++; }
        if (!hh_ok) { std::cout << " HH-MISMATCH"; fails++; }
        if (!gv_ok) { std::cout << " Giv-MISMATCH"; fails++; }
        if (qr_ok && hh_ok && gv_ok) std::cout << " ALL MATCH";
        std::cout << "\n";
    }

    std::cout << "\nFails: " << fails << "\n";
    return fails;
}
