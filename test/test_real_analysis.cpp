
/// @file test_real_analysis.cpp Test cases for real functions

#include "../src/uroboro.h"
#include <cmath>

#include "test.h"

int main(int argc, char const *argv[]) {

	std::cout << "Starting testing of uroboro library..." << std::endl;
	std::cout << "Testing real_analysis\n" << std::endl;

	std::cout.precision(8);

	test_start("umath::square(real)");

		test_tolr(umath::square(1), 1, 1);
		test_tolr(umath::square(2), 4, 2);
		test_tolr(umath::square(-1), 1, -1);
		test_tolr(umath::square(3), 9, 3);
		test_tolr(umath::square(0), 0, 0);

	test_end();


	test_start("umath::cube(real)");

		test_tolr(umath::cube(1), 1, 1);
		test_tolr(umath::cube(2), 8, 2);
		test_tolr(umath::cube(-1), -1, -1);
		test_tolr(umath::cube(3), 27, 3);
		test_tolr(umath::cube(0), 0, 0);

	test_end();


	test_start("umath::sqrt(real)");

		test_tolr_interval(umath::sqrt, std::sqrt, 0, 1, TOLERANCE, 1271351);
		test_tolr_interval(umath::sqrt, std::sqrt, 0, 1000000, TOLERANCE, 1271351);

	test_end();


	test_start("umath::cbrt(real)");

		test_tolr_interval(umath::cbrt, std::cbrt, -1, 1, TOLERANCE, 1163137);
		test_tolr_interval(umath::cbrt, std::cbrt, -10000000, 10000000, TOLERANCE, 1451319);

	test_end();


	test_start("umath::abs(real)");

		test_tolr_interval(umath::abs, std::abs, -10000000, 10000000, TOLERANCE, 1451119);

	test_end();


	test_start("umath::sgn(real)");

		test_tolr(umath::sgn(1), 1, 1);
		test_tolr(umath::sgn(2), 1, 2);
		test_tolr(umath::sgn(-1), -1, -1);
		test_tolr(umath::sgn(-3), -1, -3);
		test_tolr(umath::sgn(0), 0, 0);
		test_tolr(umath::sgn(-1.0 / 3.0), -1, -1.0 / 3.0);

	test_end();


	test_start("umath::ln(real)");

		test_tolr_interval(umath::ln, std::log, 0.00000001, 1, TOLERANCE, 1151127);
		test_tolr_interval(umath::ln, std::log, 0.00000001, 0.000001, TOLERANCE, 1151127);
		test_tolr_interval(umath::ln, std::log, 0.0001, 100000, TOLERANCE, 1151127);

	test_end();


	test_start("umath::log2(real)");

		test_tolr_interval(umath::log2, std::log2, 0.00000001, 1, TOLERANCE, 1151127);
		test_tolr_interval(umath::log2, std::log2, 0.00000001, 0.000001, TOLERANCE, 1151127);
		test_tolr_interval(umath::log2, std::log2, 0.0001, 100000, TOLERANCE, 1151127);

	test_end();


	test_start("umath::log10");

		test_tolr_interval(umath::log10, std::log10, 0.00000001, 1, TOLERANCE, 1151127);
		test_tolr_interval(umath::log10, std::log10, 0.00000001, 0.000001, TOLERANCE, 1151127);
		test_tolr_interval(umath::log10, std::log10, 0.0001, 100000, TOLERANCE, 1151127);

	test_end();


	test_start("umath::exp");

		test_tolr_interval(umath::exp, std::exp, 0, 1, TOLERANCE, 1351637);
		test_tolr_interval(umath::exp, std::exp, -10, -1, TOLERANCE, 1319673);
		test_tolr_interval(umath::exp, std::exp, 10, 20, TOLERANCE, 1137);

	test_end();


	test_start("umath::pow");

		const unsigned int N = 7;
		const unsigned int MAX_POW = 10;

		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < 100; ++j) {
				real x = 20 * i / N;
				real p = MAX_POW * j / 100;
				test_tol(umath::pow(x, p), std::pow(x, p), x, TOLERANCE, true);
			}
		}

	test_end();


	test_start("umath::powf");

		test_tolr(umath::powf(2.0, 0.5), SQRT2, 2);
		test_tolr(umath::powf(2.0, -0.5), 1.0 / SQRT2, 2);
		test_tolr(umath::powf(2, 2), 4, 2);
		test_tolr(umath::powf(3, 2), 9, 3);

	test_end();


	test_start("umath::sin");

		test_tolr_interval(umath::sin, std::sin, 0, 2 * PI, TOLERANCE, 1589167);
		test_tolr_interval(umath::sin, std::sin, 0, 100 * PI, 0.00005, 1589167);
		test_tolr_interval(umath::sin, std::sin, -100 * PI, 0, 0.00005, 1589167);

	test_end();


	test_start("umath::cos");

		test_tolr_interval(umath::cos, std::cos, 0, 2 * PI, TOLERANCE, 1589167);
		test_tolr_interval(umath::cos, std::cos, 0, 100 * PI, 0.00005, 1589167);
		test_tolr_interval(umath::cos, std::cos, -100 * PI, 0, 0.00005, 1589167);

	test_end();


	test_start("umath::tan");

		test_tolr_interval(umath::tan, std::tan, 0, PI, TOLERANCE, 1435637);

	test_end();


	test_start("umath::asin(real)");

		test_tolr_interval(umath::asin, std::asin, -0.99999, 0.99999, 0.0001);

	test_end();


	test_start("umath::acos(real)");

		test_tolr_interval(umath::acos, std::acos, -0.99999, 0.99999, 0.0001);

	test_end();


	test_start("umath::atan(real)");

		test_tolr_interval(umath::atan, std::atan, -0.5, 0.5, 0.0001, 112551);
		test_tolr_interval(umath::atan, std::atan, -1, 1, 0.0001, 12345);
		test_tolr_interval(umath::atan, std::atan, -100, 100, 0.0001, 112551);

	test_end();


	test_start("umath::sinh");

		test_tolr_interval(umath::sinh, std::sinh, 0, 1, TOLERANCE, 1726896);
		test_tolr_interval(umath::sinh, std::sinh, -5, 5, TOLERANCE, 1726896);
		test_tolr_interval(umath::sinh, std::sinh, -10, 10, TOLERANCE, 1726896);

	test_end();

	test_start("umath::cosh");

		test_tolr_interval(umath::cosh, std::cosh, 0, 1, TOLERANCE, 1726896);
		test_tolr_interval(umath::cosh, std::cosh, -5, 5, TOLERANCE, 1726896);
		test_tolr_interval(umath::cosh, std::cosh, -10, 10, TOLERANCE, 1726896);

	test_end();

	test_start("umath::tanh");

		test_tolr_interval(umath::tanh, std::tanh, 0, 1, TOLERANCE, 1726896);
		test_tolr_interval(umath::tanh, std::tanh, -5, 5, TOLERANCE, 1726896);
		test_tolr_interval(umath::tanh, std::tanh, -10, 10, TOLERANCE, 1726896);
		test_tolr_interval(umath::tanh, std::tanh, -20, 20, TOLERANCE, 1726896);

	test_end();


	test_start("umath::binomial_coeff");

		test_tol<long long, long long>(binomial_coeff(1, 1), 1, 1);
		test_tol<long long, long long>(binomial_coeff(2, 0), 1, 2);
		test_tol<long long, long long>(binomial_coeff(2, 1), 2, 2);
		test_tol<long long, long long>(binomial_coeff(3, 2), 3, 3);
		test_tol<long long, long long>(binomial_coeff(3, 1), 3, 3);
		test_tol<long long, long long>(binomial_coeff(6, 3), 20, 6);
		test_tol<long long, long long>(binomial_coeff(10, 3), 120, 10);
		test_tol<long long, long long>(binomial_coeff(16, 7), 11440, 16);
		test_tol<long long, long long>(binomial_coeff(18, 6), 18564, 18);

	test_end();


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
	else
		std::cout << "Some tests failed\n" << std::endl;


	return total_errors;
}
