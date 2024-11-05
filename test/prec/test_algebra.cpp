
/// @file test_algebra.cpp Test cases for linear algebra

#include "theoretica.h"
#include "chebyshev.h"
#include <ctime>

using namespace theoretica;
using namespace chebyshev;


// Compute the L_inf norm of any iterable structure, such as vectors or matrices.
// The norm finds the maximum element in absolute value.
template<typename Structure>
real linf_norm(const Structure& A) {

	real m = 0.0;

	for (const auto& x : A)
		m = max(m, abs(x));

	return m;
}


// Generate a random matrix with Gaussian distributed elements
mat<real> rand_mat(real m, real s, unsigned int rows, unsigned int cols) {

	mat<real> A;
	A.resize(rows, cols);

	for (auto& x : A)
		x = random::gaussian(m, s);

	return A;
}


// Generate a random symmetric matrix with Gaussian distributed elements
mat<real> rand_mat_symmetric(
	real m, real s, unsigned int rows, unsigned int cols) {

	mat<real> A = rand_mat(m, s, rows, cols);
	return (A + algebra::transpose(A)) * 0.5;
}


int main(int argc, char const *argv[]) {

	prec::settings.outputFiles = { "test/prec/test_algebra.csv" };

	prec::setup("algebra", argc, argv);

		// algebra.h

		// mat.h

		// vec.h

		// distance.h
	
	prec::terminate();
}
