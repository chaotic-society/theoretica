
///
/// @file example.cpp Example program.
/// This example may be compiled using 'make example'
///

#include "theoretica.h"
using namespace th;

#include "statistics/runstat.h"
#include "statistics/bootstrap.h"
#include <iostream>
#include <ctime>


int main() {

    PRNG g = PRNG::wyrand(time(nullptr));
    pdf_sampler gauss = pdf_sampler::gaussian(0, 1, g);
    
    vec<real> x (100);
    gauss.fill(x);
    
    std::cout << "Estimated Mean = " << stats::mean(x) << std::endl;
    std::cout << "Estimated Variance = " << stats::variance(x) << std::endl;
    std::cout << "Bootstrap = " <<
        stats::bootstrap(x, stats::mean<vec<real>>, g, 10000)
    << std::endl;
}
