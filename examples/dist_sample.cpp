
///
/// @file dist_sample.cpp Distribution sampling
///

#include "theoretica.h"

using namespace th;

#include <fstream>
#include <ctime>

int main(int argc, char const *argv[]) {

	std::cout.precision(4);

	// Ouput files
	std::ofstream f_uniform = std::ofstream("./examples/data/uniform.csv");
	std::ofstream f_gaussian = std::ofstream("./examples/data/gaussian.csv");
	std::ofstream f_exponential = std::ofstream("./examples/data/exponential.csv");
	std::ofstream f_cauchy = std::ofstream("./examples/data/cauchy.csv");
	std::ofstream f_pareto = std::ofstream("./examples/data/pareto.csv");

	// Pseudorandom Number Generator using Xoshiro256++
	PRNG g = PRNG::xoshiro(time(nullptr));

	// Discard the first 10000 values
	// to get better results
	g.discard(10000);

	// Generate 1000 values
	for (int i = 0; i < 1000; ++i) {

		// Generate a random variable following a uniform distribution
		// in the interval [0, 1]
		f_uniform << rand_real(0, 1, g) << ",\n";
		
		// Generate a random variable following a gaussian distribution
		// with mean = 0 and standard deviation = 1
		f_gaussian << rand_gaussian(0, 1, g) << ",\n";

		// Generate a random variable following an exponential distribution
		// with rate = 1
		f_exponential << rand_exponential(1, g) << ",\n";

		// Generate a random variable following a Cauchy distribution
		// with location 0 and scale 1
		f_cauchy << rand_cauchy(0, 1, g) << ",\n";
		
		// Generate a random variable following a Pareto distribution
		// with scale = 1 and shape = 2
		f_pareto << rand_pareto(1, 2, g) << ",\n";
	}

	return 0;
}
