
///
/// @file montecarlo_integral.cpp A comparison between Monte Carlo techniques
/// This example may be compiled using 'make montecarlo_integral'
///

#include <iostream>
#include "theoretica.h"
using namespace th;

#include <ctime>


int main() {


    // Initialize the pseudorandom number generator to be used
    PRNG g = PRNG::xoshiro(time(nullptr));

    // Starting number of iterations
    unsigned int N = 10;

    // The function to integrate
    real(*f)(real) = th::sin;


    // Print header
    std::cout.precision(8);
    std::cout << std::endl;
    std::cout << " N\tErr. HOM\tErr. Crude\tErr. q. HOM\tErr. q. Crude\n ";
    for (int i = 0; i < 80; ++i)
        std::cout << "-";
    std::cout << std::endl;


    // Repeat the integration for different N (exponentially growing)
    for (; N <= 100000; N *= 10) {

        // Hit-or-Miss Monte Carlo
        real hom = integral_hom(f, 0, PI/2, 1, g, N);

        // Crude Monte-Carlo
        real crude = integral_crude(f, 0, PI/2, g, N);

        // Hit-or-Miss Quasi-Monte Carlo
        real q_hom = integral_quasi_hom(f, 0, PI/2, 1, N);

        // Crude Quasi-Monte Carlo
        real q_crude = integral_quasi_crude(f, 0, PI/2, N);

        // Print absolute error
        std::cout << " "  << N << "\t"
                << th::abs(1 - hom) << "\t"
                << th::abs(1 - crude) << "\t"
                << th::abs(1 - q_hom) << "\t"
                << th::abs(1 - q_crude) << std::endl;
    }

    std::cout << std::endl;
 
    return 0;
}
