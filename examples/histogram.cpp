
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

#include <fstream>
#include <ctime>


int main() {

    // Sample size
    const unsigned int N = 1000;

    // Underlying PRNG
    PRNG g = PRNG::xoshiro(time(NULL));
            
    // Initialize a PDF sampler with a gaussian distribution
    pdf_sampler p = pdf_sampler::gaussian(0, 1, g);

    // Fill a vector with N samples
    vec_buff data(N);
    p.fill(data, N);

    // Initialize the histogram from data
    histogram h = histogram(data);

    // Print the histogram to file
    std::ofstream file("examples/histogram.dat");
    file << h;

    std::cout << "N = " << h.number() << std::endl; 
    std::cout << "Mean: " << h.mean() << std::endl; 
    std::cout << "Variance: " << h.variance() << std::endl; 
    std::cout << "Standard Deviation: " << h.stdev() << std::endl;
 
    return 0;
}
