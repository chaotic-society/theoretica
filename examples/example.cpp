
///
/// @file example.cpp Example file.
/// This example may be compiled using 'make example'
///

#include "theoretica.h"
#include <ctime>
using namespace th;


real func(real x) {
    return th::exp(-pow(x, 2)/2);
}

int main() {
    
    /*
    // Declare a 3D vector
    vec3 v = {1, 2, 3};

    // Create a 3x3 identity matrix
    mat3 A = mat3::identity();

    // Transform v by A
    vec3 w = A * v;
    */

    vec_buff samples;
    int tests = 10000;

    PRNG gen = PRNG::xoshiro(time(nullptr));
    pdf_sampler gs = pdf_sampler::gaussian(0, 1, gen);

    for(int j = 0; j < tests; j++) {
        samples.push_back(metropolis(func, gs, 1, gen));
    }

    // std::cout << "Mean: " << mean(samples) << std::endl;
    std::cout << "Variance: " << variance(samples) << std::endl;

    return 0;
}
