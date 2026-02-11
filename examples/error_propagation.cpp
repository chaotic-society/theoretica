
///
/// @file error_propagation.cpp Automatic propagation of uncertainties.
/// This example may be compiled using 'make error_propagation'
///

#include <iostream>
#include "theoretica.h"
#include <ctime>
using namespace th;


// Example function to propagate error on
template<typename NumType>
NumType f(vec<NumType> v) {

	const NumType x = v[0];
	const NumType y = v[1];
	const NumType z = v[2];

	return (x + y) * z;
}


int main() {


	// Parameters of the toy experiment
	real mu1 = 1, mu2 = 2, mu3 = 3;
	real stdev1 = 0.2, stdev2 = 0.1, stdev3 = 0.4;

	// Sample size
	const unsigned int N = 1E+06;


	// Random generators
	PRNG g = PRNG::wyrand(time(nullptr));
	pdf_sampler gauss1 = pdf_sampler::gaussian(mu1, stdev1, g);
	pdf_sampler gauss2 = pdf_sampler::gaussian(mu2, stdev2, g);
	pdf_sampler gauss3 = pdf_sampler::gaussian(mu3, stdev3, g);


	// Allocate space for 3 datasets of size N
	auto datasets = std::vector<vec<real>>(3, vec<>(N));


	// Simulate a toy experiment with Gaussian deviations
	gauss1.fill(datasets[0]);
	gauss2.fill(datasets[1]);
	gauss3.fill(datasets[2]);


	// Compute the covariance matrix
	std::cout << stats::covar_mat(datasets) << std::endl;


	std::cout << "Error:\n";

	// Propagate using the covariance matrix
	std::cout << stats::propagerr(f, datasets) << std::endl;

	// Propagate using only the standard deviation

	vec<real> means = {
		stats::mean(datasets[0]),
		stats::mean(datasets[1]),
		stats::mean(datasets[2])
	};

	vec<real> stdevs = {
		stats::stdev(datasets[0]),
		stats::stdev(datasets[1]),
		stats::stdev(datasets[2])
	};

	std::cout << stats::propagerr(f, means, stdevs) << std::endl;
}
