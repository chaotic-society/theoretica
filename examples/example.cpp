
///
/// @file example.cpp Example file.
/// This example may be compiled using 'make example'
///

#include "theoretica.h"
#include <ctime>
using namespace th;


real f(real x) {

    if(x < 0)
        return 0;

    return th::exp(-x);
}

int main() {

    vec_buff samples;
    size_t N = 1000;

    PRNG gen = PRNG::xoshiro(time(nullptr));
    pdf_sampler gauss = pdf_sampler::gaussian(0, 1, gen);

    for(unsigned int j = 0; j < N; j++)
        samples.push_back(
            metropolis(f, gauss, 2)
        );

    std::cout << "Autocorr. : " << autocorrelation(samples) << std::endl;
    std::cout << "Mean: " << mean(samples) << std::endl;
    std::cout << "Variance: " << variance(samples) << std::endl;

    return 0;
}
