#include <cmath>
#include <cstdio>
#include <iostream>
#include<limits>
#include "task_1.h"


// This function computes the machine epsilon for float, double, and long double types.
// Machine epsilon is the smallest value that, when added to 1.0, results in a value different from 1.0 due to floating-point precision limits.
// The function prints the computed values and compares them to the standard library's values.

void compute_epsilon() {
	std::cout << "--- Machine Epsilon Calculation ---\n";
	std::cout << "This program computes the machine epsilon for float, double, and long double.\n";
	std::cout << "Machine epsilon is the smallest value such that 1.0 + epsilon != 1.0.\n\n";


	float f = 1.0f;
	while ((float)(1.0f + f) != 1.0f) { f /= 2.0f; }
	f *= 2.0f;

	double d = 1.0;
	while ((double)(1.0 + d) != 1.0) { d /= 2.0; }
	d *= 2.0;

	long double l = 1.0L;
	while ((long double)(1.0L + l) != 1.0L) { l /= 2.0L; }
	l *= 2.0L;

	std::cout << "Computed machine epsilon values:\n";
	std::cout << "  float:       " << f << "\n";
	std::cout << "  double:      " << d << "\n";
	std::cout << "  long double: " << l << "\n\n";

	std::cout << "Standard library epsilon values (std::numeric_limits):\n";
	std::cout << "  float:       " << std::numeric_limits<float>::epsilon() << "\n";
	std::cout << "  double:      " << std::numeric_limits<double>::epsilon() << "\n";
	std::cout << "  long double: " << std::numeric_limits<long double>::epsilon() << "\n\n";

	std::cout << "Comparing the Values:\n";
	std::cout << "  float:       computed = " << f << ", std = " << std::numeric_limits<float>::epsilon();
	if (f == std::numeric_limits<float>::epsilon()) std::cout << " (Same!)\n";
	else std::cout << " (differ)\n";

	std::cout << "  double:      computed = " << d << ", std = " << std::numeric_limits<double>::epsilon();
	if (d == std::numeric_limits<double>::epsilon()) std::cout << " (Same!)\n";
	else std::cout << " (differ)\n";

	std::cout << "  long double: computed = " << l << ", std = " << std::numeric_limits<long double>::epsilon();
	if (l == std::numeric_limits<long double>::epsilon()) std::cout << " (Same!)\n";
	else std::cout << " (differ)\n";
    std::cout << std::endl;

    std::cout << "Theoretical values for IEEE formats:" << std::endl;
	std::cout << "  float  (2^-23):  " << std::pow(2.0f, -23) << std::endl;
	std::cout << "  double (2^-52):  " << std::pow(2.0, -52) << std::endl;
	

	std::cout << "\n--- End of Calculation ---\n";
}
