#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <functional>
#include "harm.h"

int main(int argc, char* argv[]) {
    int nthreads = 1;
    int nterms = static_cast<int>(1e8);

    for (int i = 0; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-nthreads" && i + 1 < argc) {
            nthreads = std::stoi(argv[i + 1]);
        }
        if (arg == "-nterms" && i + 1 < argc) {
            nterms = std::stoi(argv[i + 1]);
        }
    }

    std::cout << "nthreads = " << nthreads << "\n";
    std::cout << "nterms   = " << nterms << "\n";

    // Prepare data objects for each thread
    std::vector<data> params(nthreads);
    for (int i = 0; i < nthreads; i++) {
        params[i].a = 1 + nterms / nthreads * i;
        params[i].b = 1 + nterms / nthreads * (i + 1);
    }
    params[params.size() - 1].b = nterms + 1;

    // Create and launch threads
    std::vector<std::thread> threads;
    threads.reserve(nthreads);

    for (int i = 0; i < nthreads; i++) {
        threads.emplace_back(harm, std::ref(params[i]));
    }

    // Join all threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Calculate total harmonic sum
    double total = 0;
    for (auto& p : params) {
        total += p.sum;
    }

    std::cout << "Harmonic sum H(" << nterms << ") = " << total << "\n";

    return 0;
}
