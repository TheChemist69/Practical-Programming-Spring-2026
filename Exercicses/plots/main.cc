#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include "erf_func.h"
#include "sgamma_func.h"
#include "lngamma_func.h"
#include <cmath>
// Tabulated values: x, erf(x) (from standard tables)
const std::vector<std::pair<double, double>> erf_tabulated = {
    {0.0, 0.0},
    {0.1, 0.112462916},
    {0.2, 0.222702589},
    {0.3, 0.328626759},
    {0.4, 0.428392355},
    {0.5, 0.520499878},
    {0.6, 0.603856091},
    {0.7, 0.677801194},
    {0.8, 0.742100965},
    {0.9, 0.796908212},
    {1.0, 0.842700793},
    {1.5, 0.966105146},
    {2.0, 0.995322265},
    {2.5, 0.999593048},
    {3.0, 0.999977910}
};

// Tabulated values: x, gamma(x) (factorials for integer x)
const std::vector<std::pair<double, double>> gamma_tabulated = {
    {1.0, 1.0},    // 0! = 1
    {2.0, 1.0},    // 1! = 1
    {3.0, 2.0},    // 2! = 2
    {4.0, 6.0},    // 3! = 6
    {5.0, 24.0},   // 4! = 24
    {6.0, 120.0},  // 5! = 120
    {7.0, 720.0}   // 6! = 720
};

int main() {
    // Error function data
    std::ofstream erf_curve("erf_curve.dat");
    std::ofstream erf_points("erf_points.dat");
    erf_curve << std::fixed << std::setprecision(6);
    erf_points << std::fixed << std::setprecision(6);
    for (double x = -3.0; x <= 3.0; x += 0.01) {
        erf_curve << x << " " << my_erf(x) << "\n";
    }
    for (const auto& p : erf_tabulated) {
        erf_points << p.first << " " << p.second << "\n";
    }

    // Gamma function data
    std::ofstream gamma_curve("gamma_curve.dat");
    std::ofstream gamma_points("gamma_points.dat");
    gamma_curve << std::fixed << std::setprecision(6);
    gamma_points << std::fixed << std::setprecision(6);
    for (double x = -7.0; x <= 7.0; x += 0.01) {
        gamma_curve << x << " " << sgamma(x) << "\n";
    }
    for (const auto& p : gamma_tabulated) {
        gamma_points << p.first << " " << p.second << "\n";
    }

    // Log-gamma function data
    std::ofstream lngamma_curve("lngamma_curve.dat");
    std::ofstream lngamma_points("lngamma_points.dat");
    lngamma_curve << std::fixed << std::setprecision(6);
    lngamma_points << std::fixed << std::setprecision(6);
    for (double x = -5.0; x <= 5.0; x += 0.01) {
        lngamma_curve << x << " " << lngamma(x) << "\n";
    }
    for (const auto& p : gamma_tabulated) {
        lngamma_points << p.first << " " << std::log(p.second) << "\n";
    }

    std::cout << "Data written to erf_curve.dat, erf_points.dat, gamma_curve.dat, gamma_points.dat, lngamma_curve.dat, and lngamma_points.dat\n";
    return 0;
}
