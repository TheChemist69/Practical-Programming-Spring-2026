#include "input_output.h"
#include <iostream>
#include <cmath>

std::vector<double> parse_numbers(int argc, char* argv[]) {
    std::vector<double> numbers;
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-n" && i + 1 < argc) {
            numbers.push_back(std::stod(argv[i + 1]));
        }
    }
    return numbers;
}

void print_table(const std::vector<double>& numbers) {
    std::cout << "Number\tSin\t\tCos" << std::endl;
    for (auto n : numbers) {
        std::cout << n << "\t" << std::sin(n) << "\t" << std::cos(n) << std::endl;
    }
}
