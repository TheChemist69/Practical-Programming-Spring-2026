#include <iostream>
#include <iomanip>
#include "task_3.h"
#include "approx.h"

void comparing_doubles() {
    std::cout << "--- Comparing Doubles Demonstration ---\n";
    std::cout << "This demonstrates why direct comparison of floating-point numbers can fail,\n";
    std::cout << "and how to use an approximate comparison instead.\n\n";

    double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;
    double d2 = 8*0.1;

    std::cout << std::fixed << std::setprecision(17);
    std::cout << "d1 = 0.1+0.1+...+0.1 (8 times) = " << d1 << "\n";
    std::cout << "d2 = 8*0.1                 = " << d2 << "\n\n";

    std::cout << "Direct comparison (d1 == d2): " << (d1 == d2 ? "true" : "false") << "\n";
    std::cout << "Approximate comparison (approx(d1, d2)): " << (approx(d1, d2) ? "true" : "false") << "\n";

    std::cout << "\n--- End of Demonstration ---\n" << std::endl;
}
