#define UROBORO_INCLUDE_ALL
#include "../src/uroboro.h"

#include <cmath>

using namespace umath;

constexpr real TOLERANCE = 0.00000001;

std::string func_name = "";

unsigned int total_errors = 0;

unsigned int tolr_test_runs = 0;
unsigned int curr_errors = 0;
real cum_err = 0;
real max_err = 0;

void test_start(std::string f) {

	func_name = f;
	std::cout << "Testing " << func_name << " ...\n" << std::endl;

	curr_errors = 0;
	tolr_test_runs = 0;
	cum_err = 0;
	max_err = 0;
}


void test_end() {

	std::cout << "Finished testing " << func_name << std::endl;

	if(curr_errors == 0) {
		std::cout << "All tests passed successfully" << std::endl;
	} else {
		std::cout << curr_errors << " tests failed" << std::endl;
	}

	std::cout << "Absolute Error: " << cum_err << std::endl;
	std::cout << "Mean Error: " << (cum_err / (real) tolr_test_runs) << std::endl;
	std::cout << "Maximum Error: " << max_err << "\n" << std::endl;
	func_name = "";
}


template<typename T>
bool good_enough(T a, T b) {
	return umath::abs(b - a) < TOLERANCE;
}


template<typename T1, typename T2>
void test_tol(T1 evaluated, T1 expected, T2 input, bool silent = false) {

	real diff = umath::abs(evaluated - expected);
	cum_err += diff;
	max_err = umath::max(diff, max_err);
	tolr_test_runs++;

	if(!good_enough(evaluated, expected)) {

		if(!silent) {
			std::cout << "Test failed on " << func_name << ":" << std::endl;
			std::cout << "\tExpected: " << expected << std::endl;
			std::cout << "\tEvaluated: " << evaluated << std::endl;
			std::cout << "\tInput: " << input << std::endl;
			std::cout << "\tDiff: " << diff << std::endl;
		}
		
		total_errors++;
		curr_errors++;

	} else if(!silent) {
		std::cout << "Test passed with diff: " << diff << std::endl;
	}
}


void test_tolr(real evaluated, real expected, real input, bool silent = false) {
	test_tol<real, real>(evaluated, expected, input, silent);
}


void test_tolr_interval(real_function f, real_function f_exp, real a, real b, unsigned int steps = 1000) {

	std::cout << "Testing " << func_name <<
		" on interval [" << a << ", " << b << "] with " << steps << " steps" << std::endl;

	real dx = (b - a) / (real) steps;

	for (int i = 0; i <= steps; ++i) {

		real x = a + i * dx;

		test_tolr(f(x), f_exp(x), x, true);
	}
}


int main(int argc, char const *argv[]) {

	std::cout << "Starting testing of uroboro library...\n" << std::endl;

	std::cout.precision(12);


	test_start("umath::sqrt(real)");

		test_tolr(umath::sqrt(4), 2, 4);
		test_tolr(umath::sqrt(2), std::sqrt(2.0), 2);
		test_tolr(umath::sqrt(9), 3, 9);

		test_tolr_interval(umath::sqrt, std::sqrt, 0, 1);
		test_tolr_interval(umath::sqrt, std::sqrt, 0, 1000);
		test_tolr_interval(umath::sqrt, std::sqrt, 0, 10000);
		test_tolr_interval(umath::sqrt, std::sqrt, 0, 100000);
		test_tolr_interval(umath::sqrt, std::sqrt, 0, 100000);
		test_tolr_interval(umath::sqrt, std::sqrt, 0, 10000000);
		test_tolr_interval(umath::sqrt, std::sqrt, 0, 100000000);

	test_end();


	test_start("umath::ln(real)");

		test_tolr(umath::ln(E), 1, E);
		test_tolr_interval(umath::ln, std::log, 0.0001, 1000);

	test_end();


	test_start("umath::log2(real)");

		test_tolr(umath::log2(2), 1, 2);
		test_tolr(umath::log2(4), 2, 4);
		test_tolr(umath::log2(8), 3, 8);

	test_end();


	test_start("umath::log10");

		test_tolr(umath::log10(10), 1.f, 10);
		test_tolr(umath::log10(100), 2.f, 100);
		test_tolr(umath::log10(1000), 3.f, 1000);

	test_end();


	test_start("umath::exp");

		test_tolr(umath::exp(2), umath::E * umath::E, 2);
		test_tolr(umath::exp(1), umath::E, 1);

	test_end();


	test_start("umath::powf_approx");

		test_tolr(umath::powf_approx(2.f, 0.5f), umath::SQRT2, 2);

	test_end();


	test_start("umath::sin");

		// sin, cos and tan are precise
		test_tolr(umath::sin(0.5f), 0.4794255386, 0.5);
		test_tolr(umath::sin(3), 0.14112000806, 3);

	test_end();


	test_start("umath::cos");

		test_tolr(umath::cos(0.5f), 0.87758256189, 0.5);
		test_tolr(umath::cos(3), -0.9899924966, 3);

	test_end();


	test_start("umath::tan");

		test_tolr(umath::tan(0.5f), 0.54630248984, 0.5);
		test_tolr(umath::tan(3), -0.14254654307, 3);

	test_end();


	// test_start(umath::asin);

	// 	// After 0.9 gets a lot less precise
	// 	test_tolr(umath::asin(0.5), 0.5235987756, 0.5);
	// 	test_tolr(umath::asin(0.9), 1.119769515, 0.9);

	// test_end();


	// test_start(umath::acos);

	// 	// After 0.9 gets less precise
	// 	test_tolr(umath::acos(0.5), 1.0471975512, 0.5);
	// 	test_tolr(umath::acos(0.9), 0.4510268118, 0.9);

	// test_end();


	// test_start("atan(real)");

	// 	// atan is really imprecise
	// 	test_tolr(umath::atan(0.5), 0.54630248984, 0.5);
	// 	test_tolr(umath::atan(0.9), 0.78037308007, 0.9);
	// 	test_tolr_interval(umath::sqrt, std::atan, 0, 1);

	// test_end();


	test_start("umath::degrees");

		test_tolr(umath::degrees(umath::PI), 180, umath::PI);
		test_tolr(umath::degrees(umath::PI / 2.0), 90, umath::PI / 2.0);
		test_tolr(umath::degrees(umath::PI / 4.0), 45, umath::PI / 4.0);

	test_end();


	test_start("umath::radians");

		test_tolr(umath::radians(180), umath::PI, 180);
		test_tolr(umath::radians(90), umath::PI / 2.0, 90);
		test_tolr(umath::radians(45), umath::PI / 4.0, 45);

	test_end();


	// Trick test framework with fake function
	// test_start(vec3_misc);

	// 	umath::vec3 v3 = {10, 15, 20};
	// 	umath::vec4 v4 = {1, 2, 3, 4};
	// 	umath::mat4 m = umath::mat4::diagonal(2);

	// 	v4 = m * v4;

	// 	test_tolr(v4[0], 2);

	// 	v3.normalize();

	// 	test_tolr(v3.magnitude(), 1);

	// test_end();

	return total_errors;
}
