// "lsfit.h" header file.
// Ordinary least-squares fit by QR-decomposition.
#pragma once

#include <vector>
#include <functional>
#include "vector.h"
#include "matrix.h"

namespace pp {

// -----------------------------------------------------------------
// pp::lsfit -- ordinary least-squares fit of a linear combination
//   F(x) = sum_k( c_k f_k(x) )
// to measured data {x_i, y_i, dy_i}.
//
// The design matrix A has entries  A_{ik} = f_k(x_i) / dy_i ,
// and the weighted data vector  b_i = y_i / dy_i .
// The best-fit coefficients c are obtained by solving the
// overdetermined system  A c ~ b  via QR-decomposition.
//
// The covariance matrix of the fit coefficients is
//   Cov = (R^{-1})(R^{-1})^T ,
// where R is the upper-triangular factor from the QR-decomposition.
// The uncertainties are  dc_k = sqrt(Cov_{kk}) .
// -----------------------------------------------------------------

using BasisFunctions = std::vector<std::function<double(double)>>;

struct lsfit {
    BasisFunctions fs;  // basis functions used in the fit
    pp::vector c;       // best-fit coefficients
    pp::matrix cov;     // covariance matrix of the coefficients
    pp::vector dc;      // uncertainties: dc[k] = sqrt(cov[k,k])

    // Constructor: perform the least-squares fit
    explicit lsfit(
        const BasisFunctions& fs,
        const pp::vector& x,
        const pp::vector& y,
        const pp::vector& dy
    );

    // Evaluate the fit at a given point:  F(x) = sum_k( c_k f_k(x) )
    double operator()(double x) const;
};

} // namespace pp
