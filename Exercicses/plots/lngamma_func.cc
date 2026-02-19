#include "lngamma_func.h"
#include <cmath>

const double PI = 3.14159265358979323846;

double lngamma(double x) {
    if (x <= 0) return std::nan("");
    if (x < 9) return lngamma(x + 1) - std::log(x);
    return x * std::log(x + 1 / (12 * x - 1 / x / 10)) - x + std::log(2 * PI / x) / 2;
}
