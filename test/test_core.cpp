
/// @file test_core.cpp Test cases for real functions and core functionalities

#include "theoretica.h"
#include "chebyshev.h"
#include <ctime>
#include <cmath>

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {

	prec::setup("core", argc, argv);

		output::state.outputFolder = "test/";
		output::state.fieldOptions["funcName"] = output::field_options(20);
		prec::state.defaultIterations = 100'000;

		// Estimation options for real endofunctions.
		auto R_opt = prec::estimate_options<real, real>(
				prec::interval(-1E+06, 1E+06),
				prec::estimator::quadrature1D()
		);

		// Estimation options for functions defined
		// over the positive real numbers.
		auto Rplus_opt = prec::estimate_options<real, real>(
				prec::interval(0, 1E+06),
				prec::estimator::quadrature1D()
		);

		prec::estimate(
			"th::sqrt(real)",
			CAST_ENDOFUNC(th::sqrt, real),
			CAST_ENDOFUNC(std::sqrt, real),
			Rplus_opt
		);

	 	prec::estimate(
	 		"th::sqrt^2 = th::abs",
	 		[](real x) { return th::square(th::sqrt(x)); },
	 		[](real x) { return th::abs(x); },
	 		Rplus_opt
		);

		prec::estimate(
			"th::cbrt(real)",
			CAST_ENDOFUNC(th::cbrt, real),
			CAST_ENDOFUNC(std::cbrt, real),
			R_opt
		);

		prec::estimate(
			"th::cbrt^3(x) = x",
			[](real x) { return th::cube(th::cbrt(x)); },
			[](real x) { return x; },
			R_opt
		);

		prec::estimate(
			"th::root(real) (2)",
			[](real x) { return th::root(x, 2); },
			CAST_ENDOFUNC(std::sqrt, real),
			Rplus_opt
		);


		prec::estimate(
			"th::root(real) (3)",
			[](real x) { return th::root(x, 3); },
			CAST_ENDOFUNC(std::cbrt, real),
			R_opt
		);


		prec::estimate(
			"th::root(real) (4)",
			[](real x) { return th::pow(th::root(x, 4), 4); },
			[](real x) { return x; },
			Rplus_opt
		);


		prec::estimate(
			"th::isqrt(uint32_t)",
			th::isqrt<uint32_t>,
			[](real x) { return std::floor(std::sqrt(x)); },
			Rplus_opt
		);


		prec::estimate(
			"th::isqrt(uint64_t)",
			th::isqrt<uint64_t>,
			[](real x) { return std::floor(std::sqrt(x)); },
			Rplus_opt
		);


		prec::estimate(
			"th::icbrt(uint32_t)",
			th::icbrt<uint32_t>,
			[](real x) { return std::floor(std::cbrt(x)); },
			Rplus_opt
		);


		prec::estimate(
			"th::icbrt(uint64_t)",
			th::icbrt<uint64_t>,
			[](real x) { return std::floor(std::cbrt(x)); },
			Rplus_opt
		);


		prec::estimate(
			"th::abs(real)",
			CAST_ENDOFUNC(th::abs, real),
			CAST_ENDOFUNC(std::abs, real),
			R_opt
		);


		prec::estimate(
			"th::floor(real)",
			CAST_ENDOFUNC(th::floor, real),
			CAST_ENDOFUNC(std::floor, real),
			R_opt
		);


		prec::estimate(
			"th::fract(real)",
			CAST_ENDOFUNC(th::fract, real),
			[](real x) { return std::abs(std::floor(x) - x); },
			R_opt
		);


		// prec::equals("th::sgn(real)", CAST_ENDOFUNC(th::sgn), {
		// 		{1, 1},
		// 		{2, 1},
		// 		{-1, -1},
		// 		{-3, -1},
		// 		{0, 0},
		// 		{-1.0 / 3.0, -1}
		// });


		prec::estimate(
			"th::ln(real)",
			CAST_ENDOFUNC(th::ln, real),
			CAST_ENDOFUNC(std::log, real),
			prec::estimate_options<real, real>(
				prec::interval(1E-08, 1E+06),
				prec::estimator::quadrature1D()
			)
		);


		prec::estimate(
			"th::log2(real)",
			CAST_ENDOFUNC(th::log2, real),
			CAST_ENDOFUNC(std::log2, real),
			prec::estimate_options<real, real>(
				prec::interval(1E-08, 1E+06),
				prec::estimator::quadrature1D()
			)
		);


		prec::estimate(
			"th::log10(real)",
			CAST_ENDOFUNC(th::log10, real),
			CAST_ENDOFUNC(std::log10, real),
			prec::estimate_options<real, real>(
				prec::interval(1E-08, 1E+06),
				prec::estimator::quadrature1D()
			)
		);


	// 	prec::estimate(
	// 		"th::ilog2(uint32_t)",
	// 		[](real x) { return ilog2<uint32_t>(uint32_t(x)); },
	// 		[](real x) { return uint32_t(std::log2(x)); },
	// 		interval(1, MAX)
	// 	);


	// 	prec::estimate(
	// 		"th::ilog2(uint64_t)",
	// 		[](real x) { return ilog2<uint64_t>(uint64_t(x)); },
	// 		[](real x) { return uint64_t(std::log2(x)); },
	// 		interval(1, MAX));


	// 	prec::estimate(
	// 		"th::pad2(uint32_t)",
	// 		[](real x) { return pad2<uint32_t>(uint32_t(x)); },
	// 		[](real x) { return 1 << (uint32_t) ceil(std::log2(x)); },
	// 		interval(1, MAX));


	// 	prec::estimate(
	// 		"th::pad2(uint64_t)",
	// 		[](real x) { return pad2<uint64_t>(uint64_t(x)); },
	// 		[](real x) { return 1 << (uint64_t) ceil(std::log2(x)); },
	// 		interval(1, MAX));


	// 	prec::estimate(
	// 		"th::exp(real)",
	// 		CAST_ENDOFUNC(th::exp),
	// 		CAST_ENDOFUNC(std::exp),
	// 		interval(-100, 10));


	// 	prec::estimate(
	// 		"th::expm1(real)",
	// 		CAST_ENDOFUNC(th::expm1),
	// 		CAST_ENDOFUNC(std::expm1),
	// 		interval(-1, 1));


	// // test_start("th::pow");

	// // 	const unsigned int N = 7;
	// // 	const unsigned int MAX_POW = 10;

	// // 	for (unsigned int i = 0; i < N; ++i) {
	// // 		for (unsigned int j = 0; j < 100; ++j) {
	// // 			real x = 20 * i / N;
	// // 			real p = MAX_POW * j / 100;
	// // 			test_tol(th::pow(x, p), std::pow(x, p), x, TOLERANCE, true);
	// // 		}
	// // 	}

	// // test_end();

	// 	{
	// 		prec::equals("th::pow", th::pow(1, 1E+06), 1);
	// 		prec::equals("th::pow", th::pow(1, -1E+06), 1);
	// 		prec::equals("th::pow", th::pow(2, 10), (1 << 10));
	// 		prec::equals("th::pow", th::pow(10, 6), 1E+6);
	// 		prec::equals("th::pow", th::pow(E, 10) * th::pow(E, -10), 1);
	// 		prec::equals("th::pow", th::pow(1E-08, 10) * th::pow(1E-08, -10), 1);
	// 	}


	// 	{
	// 		prec::equals("th::ipow", th::ipow(1, 1E+06), 1);
	// 		prec::equals("th::ipow", th::ipow(2, 10), (1 << 10));
	// 		prec::equals("th::ipow", th::ipow(10, 6), 1E+6);
	// 	}


	// 	{
	// 		prec::equals("th::powf", th::powf(2, 0.5), th::SQRT2);
	// 		prec::equals("th::powf", th::powf(2, -0.5), 1 / th::SQRT2);
	// 		prec::equals("th::powf", th::powf(2, 2), 4);
	// 		prec::equals("th::powf", th::powf(3, 2), 9);
	// 	}


	// 	prec::estimate(
	// 		"th::sin(real)",
	// 		CAST_ENDOFUNC(th::sin),
	// 		CAST_ENDOFUNC(std::sin),
	// 		R_opt);


	// 	prec::estimate(
	// 		"th::cos(real)",
	// 		CAST_ENDOFUNC(th::cos),
	// 		CAST_ENDOFUNC(std::cos),
	// 		R_opt);


	// 	prec::estimate(
	// 		"sin^2 + cos^2 = 1",
	// 		[](real x) {
	// 			return square(th::sin(x)) + square(th::cos(x));
	// 		},
	// 		[](real x) {
	// 			return 1;
	// 		},
	// 		R_opt);


	// 	prec::estimate(
	// 		"th::tan(real)",
	// 		CAST_ENDOFUNC(th::tan),
	// 		CAST_ENDOFUNC(std::tan),
	// 		interval(-1, 1));


	// 	prec::estimate(
	// 		"th::tan(real)",
	// 		CAST_ENDOFUNC(th::tan),
	// 		CAST_ENDOFUNC(std::tan),
	// 		R_opt, 1E-6, false, 1000);

	// 	prec::equals("tan(2) = tan(2 + 100 PI)",
	// 		th::tan(2), th::tan(2 + 100 * PI));


	// 	prec::estimate(
	// 		"th::asin(real)",
	// 		CAST_ENDOFUNC(th::asin),
	// 		CAST_ENDOFUNC(std::asin),
	// 		interval(-0.999999, 0.999999),
	// 		0.0001);


	// 	prec::estimate(
	// 		"th::acos(real)",
	// 		CAST_ENDOFUNC(th::acos),
	// 		CAST_ENDOFUNC(std::acos),
	// 		interval(-0.999999, 0.999999),
	// 		0.0001);


	// 	prec::estimate(
	// 		"th::atan(real)",
	// 		CAST_ENDOFUNC(th::atan),
	// 		CAST_ENDOFUNC(std::atan),
	// 		R_opt,
	// 		0.0001);


	// 	prec::estimate(
	// 		"th::sinh(real)",
	// 		CAST_ENDOFUNC(th::sinh),
	// 		CAST_ENDOFUNC(std::sinh),
	// 		interval(-10, 10));


	// 	prec::estimate(
	// 		"th::cosh(real)",
	// 		CAST_ENDOFUNC(th::cosh),
	// 		CAST_ENDOFUNC(std::cosh),
	// 		interval(-10, 10));


	// 	prec::estimate(
	// 		"th::tanh(real)",
	// 		CAST_ENDOFUNC(th::tanh),
	// 		CAST_ENDOFUNC(std::tanh),
	// 		interval(-10, 10));

	// 	{
	// 		prec::equals("th::binomial_coeff", th::binomial_coeff(1, 1), 1, 0);
	// 		prec::equals("th::binomial_coeff", th::binomial_coeff(2, 0), 1, 0);
	// 		prec::equals("th::binomial_coeff", th::binomial_coeff(2, 1), 2, 0);
	// 		prec::equals("th::binomial_coeff", th::binomial_coeff(3, 2), 3, 0);
	// 		prec::equals("th::binomial_coeff", th::binomial_coeff(3, 1), 3, 0);
	// 		prec::equals("th::binomial_coeff", th::binomial_coeff(6, 3), 20, 0);
	// 		prec::equals("th::binomial_coeff", th::binomial_coeff(10, 3), 120, 0);
	// 		prec::equals("th::binomial_coeff", th::binomial_coeff(16, 7), 11440, 0);
	// 		prec::equals("th::binomial_coeff", th::binomial_coeff(18, 6), 18564, 0);
	// 	}


	// 	prec::equals("th::degrees(real)", CAST_ENDOFUNC(th::degrees), {
	// 		{th::PI, 180},
	// 		{th::PI/2.0, 90},
	// 		{th::PI/4.0, 45},
	// 		{0, 0},
	// 	});


	// 	prec::equals("th::radians(real)", CAST_ENDOFUNC(th::radians), {
	// 		{180, th::PI},
	// 		{90, th::PI/2.0},
	// 		{45, th::PI/4.0},
	// 		{0, 0},
	// 	});


	// 	// Square a relatively small number and check that the high bits are zero
	// 	prec::estimate("th::mul_uint128",
	// 		[](real x) {

	// 			uint64_t i = (uint64_t) x;
	// 			uint64_t r1, r2;
	// 			mul_uint128(i, i, r1, r2);

	// 			return r2;
	// 		}, [](real x) { return 0; }, interval(0, 1000));


	// 	prec::estimate("ratio::eval<real>", test_ratio, R_opt);


	// 	prec::estimate("fact<uint32_t>",
	// 		test_fact<uint32_t>, interval(1, 13));


	// 	prec::estimate("fact<uint64_t>",
	// 		test_fact<uint64_t>, interval(1, 20));


	// 	prec::estimate(
	// 		"falling_fact (0)",
	// 		[](real x) {
	// 			return falling_fact(x, 0);
	// 		},
	// 		[](real x) { return 1; },
	// 		Rplus_opt
	// 	);


	// 	prec::estimate(
	// 		"falling_fact (1)",
	// 		[](real x) {
	// 			return falling_fact(x, 1);
	// 		},
	// 		[](real x) { return x; },
	// 		Rplus_opt
	// 	);


	// 	prec::estimate(
	// 		"falling_fact (2)",
	// 		[](real x) {
	// 			return falling_fact(x, 2);
	// 		},
	// 		[](real x) { return square(x) - x; },
	// 		Rplus_opt
	// 	);


	// 	prec::estimate(
	// 		"falling_fact (3)",
	// 		[](real x) {
	// 			return falling_fact(x, 3);
	// 		},
	// 		[](real x) { return cube(x) - 3 * square(x) + 2 * x; },
	// 		Rplus_opt
	// 	);


	// 	prec::estimate(
	// 		"rising_fact (0)",
	// 		[](real x) {
	// 			return rising_fact(x, 0);
	// 		},
	// 		[](real x) { return 1; },
	// 		Rplus_opt
	// 	);


	// 	prec::estimate(
	// 		"rising_fact (1)",
	// 		[](real x) {
	// 			return rising_fact(x, 1);
	// 		},
	// 		[](real x) { return x; },
	// 		Rplus_opt
	// 	);


	// 	prec::estimate(
	// 		"rising_fact (2)",
	// 		[](real x) {
	// 			return rising_fact(x, 2);
	// 		},
	// 		[](real x) { return square(x) + x; },
	// 		Rplus_opt
	// 	);


	// 	prec::estimate(
	// 		"rising_fact (3)",
	// 		[](real x) {
	// 			return rising_fact(x, 3);
	// 		},
	// 		[](real x) { return cube(x) + 3 * square(x) + 2 * x; },
	// 		Rplus_opt
	// 	);


	// 	// Special functions


	// 	// Gamma function

	// 	// Check translation identity
	// 	prec::estimate(
	// 		"gamma (1)",
	// 		[](real x) {
	// 			return special::gamma(x);
	// 		},
	// 		[](real x) {
	// 			return special::gamma(x + 1) / x;
	// 		},
	// 		interval(0.1, 20),
	// 		10E-8, false,
	// 		prec::state.defaultIterations, prec::fail_on_rel_err
	// 	);

	// 	// Check identity with factorial
	// 	prec::estimate(
	// 		"gamma (2)",
	// 		[](real x) {
	// 			return special::gamma((real) th::floor(x));
	// 		},
	// 		[](real x) {
	// 			return (real) fact<uint64_t>((unsigned int) (th::floor(x) - 1));
	// 		},
	// 		interval(1, 20)
	// 	);


	// 	// Pi function

	// 	// Check identity with factorial
	// 	prec::estimate(
	// 		"pi (2)",
	// 		[](real x) {
	// 			return special::pi((real) th::floor(x));
	// 		},
	// 		[](real x) {
	// 			return (real) fact<uint64_t>((unsigned int) (th::floor(x)));
	// 		},
	// 		interval(1, 20)
	// 	);


	prec::terminate();
}
