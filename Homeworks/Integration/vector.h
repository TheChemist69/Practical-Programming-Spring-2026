// "vector.h" header file.
// Column-vector class for linear algebra.
#pragma once

#ifdef LONG_DOUBLE
	#define NUMBER long double
#else
	#define NUMBER double
#endif

#include <vector>
#include <initializer_list>
#include <string>
#include <functional>

namespace pp {

// -----------------------------------------------------------------
// pp::vector -- dense column-vector for numerical linear algebra
//
// Backed by std::vector<NUMBER>, where NUMBER is double by default
// (or long double if LONG_DOUBLE is defined at compile time).
//
// Supports element-wise arithmetic (+, -, *, /), scalar
// multiplication from both sides, dot product, Euclidean norm
// ||v|| = sqrt(v*v), and element-wise function application via map().
//
// Element access uses 0-based indexing:  v[i].
// Construction from initializer list enables convenient literals:
//   pp::vector v = {1.0, 2.0, 3.0};
//
// Free functions provide commutativity (scalar * vector) and a
// standalone dot() as well as approximate-equality testing.
// -----------------------------------------------------------------
class vector {
public:
    std::vector<NUMBER> data;

    // Default constructor
    vector() = default;

    // Thin wrapper around the shared Linear_Equations vector implementation.
    #pragma once

    #include "../Linear_Equations/vector.h"
    vector(std::initializer_list<NUMBER> list);
