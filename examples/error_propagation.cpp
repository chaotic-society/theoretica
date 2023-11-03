
///
/// @file error_propagation.cpp Automatic propagation of uncertainties.
/// This example may be compiled using 'make error_propagation'
///

#include <iostream>
#include "theoretica.h"
#include <ctime>
using namespace th;


// Example function to compute error on
template<typename NumType>
NumType f(vec<NumType, 3> v) {

	const NumType x = v[0];
	const NumType y = v[1];
	const NumType z = v[2];

	return (x + y) * z;
}


int main(int argc, char const *argv[]) {

	// Parameters of the toy experiment
	real mu1 = 1, mu2 = 2, mu3 = 3;
	real stdev1 = 0.2, stdev2 = 0.1, stdev3 = 0.4;

	// Random generators
	PRNG g = PRNG::xoshiro(time(nullptr));
	pdf_sampler gauss = pdf_sampler::gaussian(0, 1, g);

	std::vector<vec_buff> data;
	data.resize(3);

	// Simulate a toy experiment with Gaussian deviations
	for (int i = 0; i < 1000; ++i) {
		data[0].push_back(mu1 + gauss() * stdev1);
		data[1].push_back(mu2 + gauss() * stdev2);
		data[2].push_back(mu3 + gauss() * stdev3);
	}

	// Compute the covariance matrix
	std::cout << covar_mat<3>(data) << std::endl;

	std::cout << "Error:\n";

	// Propagate using the covariance matrix
	std::cout << propagate_err<3>(f, data) << std::endl;

	// Propagate using only the standard deviation
	std::cout << propagate_err<3>(f,
		{mean(data[0]), mean(data[1]), mean(data[2])},
		{smpl_stdev(data[0]), smpl_stdev(data[1]), smpl_stdev(data[2])}
	) << std::endl;

	return 0;
}
