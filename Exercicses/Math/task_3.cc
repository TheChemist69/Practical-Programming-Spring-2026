#include "task_2.h"
#include "task_3.h"
#include <cmath>

double lngamma(double x) {
    if (x <= 0) return NAN;
    if (x < 9) return lngamma(x + 1) - std::log(x);
    return std::log(task_2::fgamma(x));
}