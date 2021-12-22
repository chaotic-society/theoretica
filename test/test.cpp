#include "../src/uroboro.h"

#include <cmath>

using namespace umath;


// Absolute different to tolerate
constexpr real TOLERANCE = 0.00000001;

// Function name holder
std::string func_name = "";

// Total number of errors (module-wise)
unsigned int total_errors = 0;

// Number of test_tolr test runs
unsigned int tolr_test_runs = 0;

// Number of errors on current function
unsigned int curr_errors = 0;

// Cumulative error on current function
long double cum_err = 0;

// Maximum error on current function
long double max_err = 0;


// Start testing a specific function
void test_start(std::string f) {

	func_name = f;
	std::cout << "Testing " << func_name << " ...\n" << std::endl;

	curr_errors = 0;
	tolr_test_runs = 0;
	cum_err = 0;
	max_err = 0;
}


// End testing on the current function and print out information
// about the test runs
void test_end() {

	std::cout << "\nFinished testing " << func_name << std::endl;

	if(curr_errors == 0) {
		std::cout << "All tests passed successfully" << std::endl;
	} else {
		std::cout << curr_errors << " tests failed (" <<
		(curr_errors / (real) tolr_test_runs * 100) << "%)" << std::endl;
	}

	std::cout << "Mean Error: " << (cum_err / (real) tolr_test_runs) << std::endl;
	std::cout << "Maximum Error: " << max_err << "\n\n" << std::endl;
	func_name = "";
}


// Check whether the given values differ only by a tolerance value or less
template<typename T>
bool good_enough(T a, T b, real tolerance = TOLERANCE) {
	return umath::abs(b - a) < tolerance;
}


// Compare a function to an expected value
template<typename T1, typename T2>
void test_tol(T1 evaluated, T1 expected, T2 input, real tolerance = TOLERANCE, bool silent = false) {

	real diff = umath::abs(evaluated - expected);
	cum_err += diff;
	max_err = umath::max(diff, max_err);
	tolr_test_runs++;

	if(!good_enough(evaluated, expected, tolerance)) {

		if(!silent) {
			std::cout << "\tTest failed on " << func_name << ":" << std::endl;
			std::cout << "\t\tExpected: " << expected << std::endl;
			std::cout << "\t\tEvaluated: " << evaluated << std::endl;
			std::cout << "\t\tInput: " << input << std::endl;
			std::cout << "\t\tDiff: " << diff << std::endl;
		}
		
		total_errors++;
		curr_errors++;

	} else if(!silent) {
		std::cout << "\tTest passed with diff: " << diff << std::endl;
	}
}


void test_tolr(real evaluated, real expected, real input, real tolerance = TOLERANCE, bool silent = false) {
	test_tol<real, real>(evaluated, expected, input, tolerance, silent);
}



// Test a real function on an interval
void test_tolr_interval(real_function f, real_function f_exp, real a, real b,
	real tolerance = TOLERANCE, unsigned int steps = 1000) {

	std::cout << "\tTesting on interval [" << a << ", " << b << "] with " << steps << " steps" << std::endl;

	real dx = (b - a) / (real) steps;

	for (int i = 0; i <= steps; ++i) {

		real x = a + i * dx;

		test_tolr(f(x), f_exp(x), x, tolerance, true);
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
		test_tolr_interval(umath::ln, std::log, 0.00000001, 1);
		test_tolr_interval(umath::ln, std::log, 0.00000001, 0.000001);
		test_tolr_interval(umath::ln, std::log, 0.0001, 1000);

	test_end();


	test_start("umath::log2(real)");

		test_tolr(umath::log2(2), 1, 2);
		test_tolr(umath::log2(4), 2, 4);
		test_tolr(umath::log2(8), 3, 8);
		test_tolr_interval(umath::log2, std::log2, 0.00000001, 1);
		test_tolr_interval(umath::log2, std::log2, 0.0001, 1000);

	test_end();


	test_start("umath::log10");

		test_tolr(umath::log10(10), 1.f, 10);
		test_tolr(umath::log10(100), 2.f, 100);
		test_tolr(umath::log10(1000), 3.f, 1000);
		test_tolr_interval(umath::log10, std::log10, 0.00000001, 1);
		test_tolr_interval(umath::log10, std::log10, 0.0001, 1000);

	test_end();


	test_start("umath::exp");

		test_tolr(umath::exp(2), E * E, 2);
		test_tolr(umath::exp(1), E, 1);
		// test_tolr_interval(umath::exp, std::exp, 0.00000001, 1);
		// test_tolr_interval(umath::exp, std::exp, 1, 10);

	test_end();


	test_start("umath::powf_approx");

		test_tolr(umath::powf_approx(2.0, 0.5), SQRT2, 2);
		test_tolr(umath::powf_approx(2, 2), 4, 2);
		test_tolr(umath::powf_approx(3, 2), 9, 3);

	test_end();


	test_start("umath::sin");

		test_tolr(umath::sin(0.5f), 0.4794255386, 0.5);
		test_tolr(umath::sin(3), 0.14112000806, 3);
		test_tolr_interval(umath::sin, std::sin, 0, 2 * PI, TOLERANCE, 100000);
		test_tolr_interval(umath::sin, std::sin, 0, 10 * PI);
		test_tolr_interval(umath::sin, std::sin, -10 * PI, 0);
		test_tolr_interval(umath::sin, std::sin, 0, 100 * PI);
		test_tolr_interval(umath::sin, std::sin, -100 * PI, 0);
		test_tolr_interval(umath::sin, std::sin, 0, 1000 * PI);
		test_tolr_interval(umath::sin, std::sin, -1000 * PI, 0);

	test_end();


	test_start("umath::cos");

		test_tolr(umath::cos(0.5f), 0.87758256189, 0.5);
		test_tolr(umath::cos(3), -0.9899924966, 3);
		test_tolr_interval(umath::cos, std::cos, 0, 2 * PI, TOLERANCE, 100000);
		test_tolr_interval(umath::cos, std::cos, 0, 10 * PI);
		test_tolr_interval(umath::cos, std::cos, -10 * PI, 0);
		test_tolr_interval(umath::cos, std::cos, 0, 100 * PI);
		test_tolr_interval(umath::cos, std::cos, -100 * PI, 0);
		test_tolr_interval(umath::cos, std::cos, 0, 1000 * PI);
		test_tolr_interval(umath::cos, std::cos, -1000 * PI, 0);

	test_end();


	test_start("umath::tan");

		test_tolr(umath::tan(0.5f), 0.54630248984, 0.5);
		test_tolr(umath::tan(3), -0.14254654307, 3);
		test_tolr_interval(umath::tan, std::tan, 0, PI, TOLERANCE, 300000);

	test_end();


	// test_start("umath::asin(real)");

	// 	test_tolr_interval(umath::asin, std::asin, 0, 1, 0.0001);
	// 	test_tolr_interval(umath::asin, std::asin, -1, 0, 0.0001);

	// test_end();


	// test_start("umath::acos(real)");

	// 	test_tolr_interval(umath::acos, std::acos, 0, 1, 0.0001);
	// 	test_tolr_interval(umath::acos, std::acos, -1, 0, 0.0001);

	// test_end();


	// test_start("umath::atan(real)");

	// 	// Allow greater tolerance on atan
	// 	test_tolr(umath::atan(0.5), std::atan(0.5), 0.5, 0.0001);
	// 	test_tolr(umath::atan(0.9), std::atan(0.9), 0.9, 0.0001);
	// 	test_tolr_interval(umath::atan, std::atan, 0, 1, 0.0001);
	// 	test_tolr_interval(umath::atan, std::atan, -1, 0, 0.0001);

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

	if(total_errors == 0)
		std::cout << "All tests on all functions and modules successfully passed\n" << std::endl;


	return total_errors;
}
