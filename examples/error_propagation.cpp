
///
/// @file error_propagation.cpp Automatic error propagation using automatic differentiation
///

#include "theoretica.h"
using namespace th;

#include <iostream>


// Direct sum of the errors
template<unsigned int N>
real product_sum(vec<N> g, vec<N> v) {

	real res = 0;
	for (unsigned int i = 0; i < N; ++i)
		res += th::abs(g.at(i) * v.at(i));

	return res;
}

// Sum the errors under quadrature
template<unsigned int N>
real product_sum_quad(vec<N> g, vec<N> v) {

	real res = 0;
	for (unsigned int i = 0; i < N; ++i)
		res += square(g.at(i) * v.at(i));

	return th::sqrt(res);
}



// Example function to compute error on
template<typename NumType>
NumType theta(vec<2, NumType> v) {

	const NumType r = v[0];
	const NumType d = v[1];

	return th::atan(r / (d - 1));
}



int main(int argc, char const *argv[]) {

	real r, d;
	real delta_r, delta_d;

	// Read data and errors from stdin
	std::cout << "r = ";
	std::cin >> r;

	std::cout << "delta_r = ";
	std::cin >> delta_r;

	std::cout << "d = ";
	std::cin >> d;

	std::cout << "delta_d = ";
	std::cin >> delta_d;

	// Construct data and error as vectors
	vec2 data = {r, d};
	vec2 err = {delta_r, delta_d};

	// Compute the partial derivatives of the function
	vec2 gradient = th::gradient(theta, data);

	// Propagate the error on the function by adding
	// the absolute values of the partial derivative
	// times the error on the variable
	// (either in direct sum or quadrature)

	std::cout << "Total error (direct sum): " <<
		product_sum(gradient, err) << std::endl;

	std::cout << "Total error (sum under quadrature): " <<
		product_sum_quad(gradient, err) << std::endl;

	return 0;
}
