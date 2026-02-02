
#include <iostream>
#include <iomanip>
#include <cmath>
#include "task_2.h"

// This function demonstrates how tiny values (smaller than machine epsilon) affect floating-point addition order.
void tiny_epsilon() {
	double epsilon = std::pow(2, -52); // IEEE double precision epsilon
	double tiny = epsilon / 2;
	double a = 1 + tiny + tiny;
	double b = tiny + tiny + 1;

	std::cout << "--- Tiny Epsilon Demonstration ---\n";
	std::cout << "This program shows how adding very small values (less than machine epsilon) to 1.0 can behave differently depending on the order of addition.\n";
	std::cout << "We use epsilon = 2^-52 (double precision) and tiny = epsilon/2.\n\n";

	std::cout << "Values used:\n";
	std::cout << "  epsilon (2^-52):      " << std::setprecision(17) << epsilon << "\n";
	std::cout << "  tiny (epsilon/2):     " << tiny << "\n\n";

	std::cout << "Computed results:\n";
	std::cout << std::fixed << std::setprecision(17);
	std::cout << "a=  1 + tiny + tiny =     " << a << "\n";
	std::cout << "b=  tiny + tiny + 1 =     " << b << "\n\n";

	std::cout << "Comparisons (using above definitions):\n";
	std::cout << "  a == b ?              " << (a == b ? "true" : "false") << "\n";
	std::cout << "  a > 1  ?              " << (a > 1  ? "true" : "false") << "\n";
	std::cout << "  b > 1  ?              " << (b > 1  ? "true" : "false") << "\n";

	std::cout << "\nThis demonstrates that floating-point addition is not always associative when tiny values are involved.\n";
	std::cout << "--- End of Demonstration ---\n";
}