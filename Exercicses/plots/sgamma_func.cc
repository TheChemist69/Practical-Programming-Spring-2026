#include "sgamma_func.h"
#include <cmath>

const double PI = 3.14159265358979323846;

double sgamma(double x) {
    if (x < 0) return PI / std::sin(PI * x) / sgamma(1 - x);
    if (x < 9) return sgamma(x + 1) / x;
    double lnsgamma = std::log(2 * PI) / 2 + (x - 0.5) * std::log(x) - x
        + (1.0 / 12) / x - (1.0 / 360) / (x * x * x) + (1.0 / 1260) / (x * x * x * x * x);
    return std::exp(lnsgamma);
}
