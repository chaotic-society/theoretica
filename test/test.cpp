
#include "../src/uroboro.h"
#include <cmath>
#include "test.h"


int main(int argc, char const *argv[]) {

	std::cout << "Starting testing of uroboro library...\n" << std::endl;

	std::cout.precision(12);


	test_start("umath::sqrt(real)");

		test_tolr(umath::sqrt(4), 2, 4);
		test_tolr(umath::sqrt(2), std::sqrt(2.0), 2);
		test_tolr(umath::sqrt(9), 3, 9);

		test_tolr_interval(umath::sqrt, std::sqrt, 0, 1, TOLERANCE, 100000);
		test_tolr_interval(umath::sqrt, std::sqrt, 0, 1000);
		test_tolr_interval(umath::sqrt, std::sqrt, 0, 10000);
		test_tolr_interval(umath::sqrt, std::sqrt, 0, 100000);
		test_tolr_interval(umath::sqrt, std::sqrt, 0, 10000000);
		test_tolr_interval(umath::sqrt, std::sqrt, 0, 100000000);

	test_end();


	test_start("umath::cbrt(real)");

		test_tolr(umath::cbrt(8), 2, 8);
		test_tolr(umath::cbrt(2), std::cbrt(2.0), 2);
		test_tolr(umath::cbrt(27), 3, 27);

		test_tolr_interval(umath::cbrt, std::cbrt, 0, 1, TOLERANCE, 100000);
		test_tolr_interval(umath::cbrt, std::cbrt, 0, 1000);
		test_tolr_interval(umath::cbrt, std::cbrt, 0, 10000);
		test_tolr_interval(umath::cbrt, std::cbrt, 0, 100000);
		test_tolr_interval(umath::cbrt, std::cbrt, 0, 10000000);
		test_tolr_interval(umath::cbrt, std::cbrt, 0, 100000000);

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
		test_tolr_interval(umath::exp, std::exp, 0.00000001, 1);

	test_end();


	test_start("umath::powf_approx");

		test_tolr(umath::powf_approx(2.0, 0.5), SQRT2, 2);
		test_tolr(umath::powf_approx(2, 2), 4, 2);
		test_tolr(umath::powf_approx(3, 2), 9, 3);

	test_end();


	test_start("umath::sin");

		test_tolr(umath::sin(0.5f), 0.4794255386, 0.5);
		test_tolr(umath::sin(3), 0.14112000806, 3);
		test_tolr_interval(umath::sin, std::sin, 0, PI, TOLERANCE, 100000);
		// test_tolr_interval(umath::sin, std::sin, 0, 10 * PI);
		// test_tolr_interval(umath::sin, std::sin, -10 * PI, 0);
		// test_tolr_interval(umath::sin, std::sin, 0, 100 * PI);
		// test_tolr_interval(umath::sin, std::sin, -100 * PI, 0);
		// test_tolr_interval(umath::sin, std::sin, 0, 1000 * PI);
		// test_tolr_interval(umath::sin, std::sin, -1000 * PI, 0);

	test_end();


	test_start("umath::cos");

		test_tolr(umath::cos(0.5f), 0.87758256189, 0.5);
		test_tolr(umath::cos(3), -0.9899924966, 3);
		test_tolr_interval(umath::cos, std::cos, 0, PI2, TOLERANCE, 100000);
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
		test_tolr_interval(umath::tan, std::tan, 0, PI, TOLERANCE, 333333);

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
	// 	test_tolr(umath::atan(0.5), std::atan(0.5), 0.5);
	// 	test_tolr(umath::atan(0.9), std::atan(0.9), 0.9);
	// 	test_tolr_interval(umath::atan, std::atan, 0, 1);
	// 	test_tolr_interval(umath::atan, std::atan, -1, 0);

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


	if(total_errors == 0)
		std::cout << "All tests on all functions and modules successfully passed\n" << std::endl;
	else
		std::cout << "Some tests failed\n" << std::endl;


	return total_errors;
}
