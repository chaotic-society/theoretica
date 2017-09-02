#include <iostream>
#define UROBORO_DOUBLE_PRECISION
#include "uroboro.h"
#include "utility.h"

using namespace uroboro;

constexpr real TOLERANCE = 0.01;

bool good_enough(real a, real b) {
	return abs(b - a) < TOLERANCE;
}

real square(real x) {
	return x * x;
}

void print_real(real x) {
	std::cout << x << std::endl;
}

void test(real a, real b) {
	if(!good_enough(a, b)) {
		std::cout << "Test not passed" << std::endl;
		std::cout << "\tResult " << a << std::endl;
		std::cout << "\tExpected " << b << std::endl;
	} else
		std::cout << "Test passed" << std::endl;
}


int main(int argc, char const *argv[]) {

	set_output_prec(8);
	test(sqrt(4), 2);
	test(sqrt(2), SQRT2);

	test(log(E), 1);
	test(log(E * E), 2);

	test(log2(2), 1);
	test(log2(4), 2);
	test(log2(8), 3);

	// log10 is safe only for x in range 0-20
	test(log10(10), 1.f);
	test(log10(100), 2.f);
	test(log10(1000), 3.f);

	test(exp(2), E * E);
	test(exp(1), E);

	// powf is safe only for x in range 0-1
	test(powf(1.f, 0.5f), 1.f);

	// sin, cos and tan are precise
	test(sin(0.5f), 0.4794255386);
	test(sin(3), 0.14112000806);

	test(cos(0.5f), 0.87758256189);
	test(cos(3), -0.9899924966);

	test(tan(0.5f), 0.54630248984);
	test(tan(3), -0.14254654307);

	// After 0.9 gets less precise
	test(asin(0.5f), 0.5235987756);
	test(asin(0.9f), 1.119769515);

	// After 0.9 gets less precise
	test(acos(0.5f), 1.0471975512);
	test(acos(0.9f), 0.4510268118);

	// atan is really imprecise
	test(atan(0.5f), 0.54630248984);
	test(atan(0.9f), 0.78037308007);

	test(degrees(50), 2864.789f);
	test(radians(50), 0.8726646f);

	vec4();

	return 0;
}
