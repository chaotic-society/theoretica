
///
/// @file histogram.cpp Histogram usage example.
/// This example may be compiled using 'make histogram'.
///
/// The resulting histogram file can be easily be visualized
/// using the gnuplot command:
/// plot "examples/histogram.dat" with boxes
///

#include "theoretica.h"
using namespace th;

#include <ctime>
#include <fstream>
#include <iostream>


int main() {

    // Sample size
    const unsigned int N = 1E+06;

    // Create a PRNG with a seed
    PRNG g = PRNG::xoshiro(time(nullptr));
            
    // Initialize a p.d.f. sampler with a Gaussian distribution
    pdf_sampler p = pdf_sampler::gaussian(0, 1, g);

    // Fill a vector with N samples
    vec<real> data(N);
    p.fill(data);

    // Initialize the histogram from data
    histogram h = histogram(data);

    // If no additional options are passed,
    // the histogram class uses default parameters
    // which work for most use cases.

    // Print the histogram to file
    std::ofstream file("examples/histogram.dat");
    file << h;

    // Print fundamental statistics of the sample
    std::cout << "N = " << h.number() << std::endl; 
    std::cout << "Mean: " << stats::mean(h) << std::endl; 
    std::cout << "Variance: " << stats::variance(h) << std::endl; 
    std::cout << "Standard Deviation: " << stats::stdev(h) << std::endl;
}
