
/// @file test_algebra.cpp Test cases for linear algebra

#include "theoretica.h"
#include "chebyshev/prec.h"
#include <ctime>

using namespace theoretica;
using namespace chebyshev;


// Generate a random matrix with elements inside the interval [a, b]
template<unsigned int N, unsigned int M>
mat<N, M> rand_mat(real a, real b, PRNG& g) {

	mat<N, M> A;
	for (unsigned int i = 0; i < N; ++i)
		for (unsigned int j = 0; j < M; ++j)
			A.iat(i, j) = rand_uniform(a, b, g);

	return A;
}


// Test mat<N, N>::inverse()
template<unsigned int N>
prec::estimate_result test_matrix_inverse(interval k, Real tol, unsigned int n) {

	Real max = 0;
	Real sum = 0;
	Real sum2 = 0;

	PRNG g = PRNG::xoshiro(time(nullptr));
	g.discard(1000);

	for (unsigned int i = 0; i < n; ++i) {

		mat<N, N> A = rand_mat<N, N>(k.a, k.b, g);

		// Skip singular matrices
		if(th::abs(A.det()) <= MACH_EPSILON) {
			if(i) i--;
			continue;
		}

		// Resulting matrix expected to be identity
		mat<N, N> R = A * A.inverse();

		for (unsigned int j = 0; j < N; ++j) {
			for (unsigned int k = 0; k < N; ++k) {

				Real diff = th::abs(R.iat(j, k) - kronecker_delta(j, k));

				sum += diff;
				sum2 += square(diff);
				if(diff > max)
					max = diff;
			}
		}

	}

	prec::estimate_result res;
	res.max_err = max;
	res.abs_err = sum / n;
	res.rms_err = th::sqrt(sum2) / n;
	res.mean_err = sum / N / n;

	// Undefined relative error
	res.rel_err = 0;

	if(res.max_err > tol)
		res.failed = true;

	return res;
}



int main(int argc, char const *argv[]) {

	prec::state.outputFolder = "test/";

	std::vector<interval> intervals = {
		interval(-1, 1),
		interval(-10000000, 10000000)
	};

	prec::setup("algebra");

		prec::estimate("mat3::inverse", test_matrix_inverse<3>, intervals);
		prec::estimate("mat4::inverse", test_matrix_inverse<4>, intervals);
		prec::estimate("mat10::inverse", test_matrix_inverse<10>, intervals);

	prec::terminate();
}
