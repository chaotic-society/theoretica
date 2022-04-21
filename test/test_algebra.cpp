
/// @file test_algebra.cpp Test cases for linear algebra

#include "../src/uroboro.h"
#include <cmath>

#include "test.h"
#include <ctime>


int main(int argc, char const *argv[]) {

	std::cout << "Starting testing of uroboro library..." << std::endl;
	std::cout << "Testing algebra/mat\n" << std::endl;

	std::cout.precision(12);


	test_start("umath::mat::inverse");

		unsigned int N = 1000000;

		std::cout << "\tTesting on " << N << " random matrices" << std::endl;

		PRNG g = PRNG::linear_congruential(time(nullptr));

		for (int i = 0; i < N; ++i) {
			
			mat4 A;

			// Generate a random matrix
			for (int j = 0; j < 4; ++j) {
				for (int k = 0; k < 4; ++k) {
					A.iat(j, k) = rand_real(-1000000, 1000000, g);
				}
			}

			// Skip singular matrices
			if(A.det() == 0) {
				i--;
				continue;
			}

			mat4 Ainv = A.inverse();
			mat4 res = A * Ainv;

			// Check that all entries are zero
			for (int j = 0; j < 4; ++j) {
				for (int k = 0; k < 4; ++k) {
					test_tolr(
						res.iat(j, k),
						kronecker_delta(j, k),
						A.iat(j, k),
						TOLERANCE, true);
				}
			}

		}

	test_end();

	if(total_errors == 0)
		std::cout << "All tests on all functions and modules successfully passed\n" << std::endl;
	else
		std::cout << "Some tests failed\n" << std::endl;


	return total_errors;
}
