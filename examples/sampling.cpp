
///
/// @file sampling.cpp Distribution sampling example.
/// This example may be compiled using 'make sampling'
///

#include "theoretica.h"
#include <iostream>
#include <fstream>
#include <ctime>

using namespace th;


int main() {

	// Sample size
	const unsigned int N = 10000;

	// Output files
	std::ofstream f_uniform = std::ofstream("./examples/uniform.csv");
	std::ofstream f_gaussian = std::ofstream("./examples/gaussian.csv");
	std::ofstream f_exponential = std::ofstream("./examples/exponential.csv");
	std::ofstream f_cauchy = std::ofstream("./examples/cauchy.csv");

	// Pseudorandom Number Generator using Xoshiro256++
	random::XoshiroPrng g (time(nullptr));

	// Generate N values
	for (unsigned int i = 0; i < N; ++i) {

		// Generate a random variable following a uniform distribution
		// in the interval [0, 1]
		f_uniform << random::uniform(0, 1, g) << ",\n";
		
		// Generate a random variable following a Gaussian distribution
		// with mean = 0 and standard deviation = 1
		f_gaussian << random::gaussian(0, 1, g) << ",\n";

		// Generate a random variable following an exponential distribution
		// with rate = 1
		f_exponential << random::exponential(1, g) << ",\n";

		// Generate a random variable following a Cauchy distribution
		// with location 0 and scale 1
		f_cauchy << random::cauchy(0, 1, g) << ",\n";
	}
}
