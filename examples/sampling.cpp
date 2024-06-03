
///
/// @file sampling.cpp Distribution sampling example.
/// This example may be compiled using 'make sampling'
///

#include "theoretica.h"
#include <iostream>
#include <fstream>
#include <ctime>

using namespace th;


int main(int argc, char const *argv[]) {

	// Sample size
	const unsigned int N = 10000;

	// Output files
	std::ofstream f_uniform = std::ofstream("./examples/data/uniform.csv");
	std::ofstream f_gaussian = std::ofstream("./examples/data/gaussian.csv");
	std::ofstream f_exponential = std::ofstream("./examples/data/exponential.csv");
	std::ofstream f_cauchy = std::ofstream("./examples/data/cauchy.csv");

	// Pseudorandom Number Generator using Xoshiro256++
	PRNG g = PRNG::xoshiro(time(nullptr));

	// Discard the first 10000 values
	// to get better results from the PRNG
	g.discard(10000);

	// Generate N values
	for (unsigned int i = 0; i < N; ++i) {

		// Generate a random variable following a uniform distribution
		// in the interval [0, 1]
		f_uniform << rand_uniform(0, 1, g) << ",\n";
		
		// Generate a random variable following a Gaussian distribution
		// with mean = 0 and standard deviation = 1
		f_gaussian << rand_gaussian(0, 1, g) << ",\n";

		// Generate a random variable following an exponential distribution
		// with rate = 1
		f_exponential << rand_exponential(1, g) << ",\n";

		// Generate a random variable following a Cauchy distribution
		// with location 0 and scale 1
		f_cauchy << rand_cauchy(0, 1, g) << ",\n";
	}
}
