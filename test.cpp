#include <iostream>
#define UROBORO_DOUBLE_PRECISION
#include "uroboro.h"
#include "utility.h"

using namespace uroboro;

constexpr real TOLERANCE = 0.001;

bool good_enough(real a, real b) {
	return abs(b - a) < TOLERANCE;
}

real square(real x) {
	return x * x;
}

void print_real(real x) {
	std::cout << x << std::endl;
}


int main(int argc, char const *argv[]) {

	set_output_prec(8);

	return 0;
}
