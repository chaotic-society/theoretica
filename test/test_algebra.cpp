
#include "../src/uroboro.h"
#include <cmath>

#include "test.h"


int main(int argc, char const *argv[]) {

	std::cout << "Starting testing of uroboro library..." << std::endl;
	std::cout << "Testing algebra/mat\n" << std::endl;

	std::cout.precision(12);


	test_start("umath::mat::inverse");

		mat4 A4 = mat4(3);
		A4.at(0, 1) = 2;
		A4.at(1, 0) = 1;
		A4.at(3, 1) = 5;

		// test_equal(A4.inverse() * A4, mat4::identity(), A4);

	test_end();

	if(total_errors == 0)
		std::cout << "All tests on all functions and modules successfully passed\n" << std::endl;
	else
		std::cout << "Some tests failed\n" << std::endl;


	return total_errors;
}
