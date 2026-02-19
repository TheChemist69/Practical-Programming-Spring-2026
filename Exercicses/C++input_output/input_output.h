#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include <vector>
#include <string>

std::vector<double> parse_numbers(int argc, char* argv[]);
void print_table(const std::vector<double>& numbers);

#endif // INPUT_OUTPUT_H
