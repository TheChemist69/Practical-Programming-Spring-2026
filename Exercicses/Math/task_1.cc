#include <iostream>
#include <cmath>
#include "task_1.h"

void task_1( void) {
    double sqrt2 = std::sqrt(2.0);                              
    std::cout << "The square root of 2 is: " << sqrt2 << "\n" ; 

    double power = std::pow(3, 1.0/5.0);
    std::cout << "3 raised to the power of 1/5 is: " << power << "\n" ;

    double exp = std::exp(M_PI);
    std::cout << "e raised to the power of pi is: " << exp << "\n" ;

    double Pi_e = std::pow(M_PI, std::exp(1.0));
    std::cout << "pi raised to the power of e is: " << Pi_e << "\n" ;
}