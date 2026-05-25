#include "harm.h"

void harm(data& arg) {
    int a = arg.a;
    int b = arg.b;
    double sum = 0;
    for (int i = a; i < b; i++) {
        sum += 1.0 / i;
    }
    arg.sum = sum;
}
