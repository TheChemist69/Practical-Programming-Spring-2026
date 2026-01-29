#include <iostream>
#include "task_1.h"
#include "task_2.h"
#include "task_3.h"

int main() {
     // task 1, doing simple math operations
     std::cout << "task 1: Simple Math Operations" << "\n";
     task_1();

     // printing a separator
     std::cout << "------------------------" << "\n";

     // task 2, computing the gamma function
     std::cout << "task 2: Computing the Gamma Function" << "\n";
     double test_values[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 500.0, -2.0};

     for (double x : test_values) {
          double result = task_2::fgamma(x);
          std::cout << "Gamma(" << x << ") = " << result << "\n";
     }

     // printing a separator
     std::cout << "------------------------" << "\n";
     // task 3, computing the natural logarithm of the gamma function
     std::cout << "task 3: Computing the Natural Logarithm of the Gamma Function" << "\n";
     for (double x : test_values) {
          double result = lngamma(x);
          std::cout << "ln(Gamma(" << x << ")) = " << result << "\n";
     }




     return 0;
}