#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "input_output.h"
#include <cstdlib>
#include <iterator>

int main(int argc, char* argv[]) {
    std::string infile = "", outfile = "";
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--input" && i + 1 < argc) infile = argv[i + 1];
        if (arg == "--output" && i + 1 < argc) outfile = argv[i + 1];
    }
    if (!infile.empty() && !outfile.empty()) {
        std::ifstream myinput(infile);
        std::ofstream myoutput(outfile);
        double x;
        if (myinput.is_open() && myoutput.is_open()) {
            myoutput << "Number\tSin\t\tCos" << std::endl;
            while (myinput >> x) {
                myoutput << x << "\t" << std::sin(x) << "\t" << std::cos(x) << std::endl;
            }
        } else {
            std::cerr << "Error opening files: " << infile << " " << outfile << std::endl;
            return EXIT_FAILURE;
        }
        myinput.close();
        myoutput.close();
        return EXIT_SUCCESS;
    }
    // Fallback to previous behavior
    std::vector<double> numbers = parse_numbers(argc, argv);
    if (numbers.empty()) {
        double x;
        while (std::cin >> x) {
            numbers.push_back(x);
        }
    }
    print_table(numbers);
    return EXIT_SUCCESS;
}