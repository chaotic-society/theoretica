
/// @file test_algebra.cpp Test cases for linear algebra

#include "../src/theoretica.h"
#include <cmath>

#include "test.h"
#include <ctime>


constexpr real FAIL_THRESHOLD_PERC = 1E-04;


int main(int argc, char const *argv[]) {

	std::cout << "Starting testing of theoretica library..." << std::endl;
	std::cout << "Testing algebra/mat\n" << std::endl;

	std::cout.precision(12);


	test_start("th::mat::inverse");

		unsigned int N = 1000000;

		std::cout << "\tTesting on " << N << " random matrices" << std::endl;

		PRNG g = PRNG::xoshiro(time(nullptr));
		g.discard(1000);

		for (unsigned int i = 0; i < N; ++i) {
			
			mat4 A;

			// Generate a random matrix
			for (unsigned int j = 0; j < 4; ++j) {
				for (int k = 0; k < 4; ++k) {
					A.iat(j, k) = rand_uniform(-1000000, 1000000, g);
				}
			}

			// Skip singular matrices
			if(A.det() == 0) {
				i--;
				continue;
			}

			mat4 Ainv = A.inverse();
			mat4 res = A * Ainv;

			// Check that all entries are zero except on the diagonal
			for (unsigned int j = 0; j < 4; ++j) {
				for (unsigned int k = 0; k < 4; ++k) {
					test_tolr(
						res.iat(j, k),
						kronecker_delta(j, k),
						A.iat(j, k),
						TOLERANCE, true);
				}
			}

		}

		if(curr_errors / (real) tolr_test_runs <= FAIL_THRESHOLD_PERC)
			total_errors -= curr_errors;

	test_end();

	if(total_errors == 0)
		std::cout << "All tests on all functions and modules successfully passed\n" << std::endl;
	else
		std::cout << "Some tests failed\n" << std::endl;


	return total_errors;
}
