
/// @file test_core.cpp Test cases for real functions and core functionalities

#include "theoretica.h"
#include "chebyshev/prec.h"
#include <ctime>
#include <cmath>

using namespace chebyshev;
using namespace theoretica;

prec::estimate_result test_ratio(interval k, Real tol, unsigned int n) {

	Real max = 0;
	Real sum = 0;
	Real sum2 = 0;

	PRNG g = PRNG::xoshiro(time(nullptr));
	g.discard(1000);

	for (unsigned int i = 0; i < n; ++i) {

		real a, b, c, d;
		a = rand_uniform(k.a, k.b, g);
		c = rand_uniform(k.a, k.b, g);
		b = rand_uniform(k.a, k.b, g);
		d = rand_uniform(k.a, k.b, g);

		ratio<real, real> AB = ratio<real, real>(a, b);
		ratio<real, real> CD = ratio<real, real>(c, d);

		real diff = th::abs(((a * c) / (b * d)) - real(AB * CD));
		
		if(max < diff)
			max = diff;

		sum += diff;
		sum2 += square(diff);
	}

	prec::estimate_result p;
	p.max_err = max;
	p.abs_err = sum / n;
	p.rms_err = th::sqrt(sum2) / n;
	p.mean_err = sum / n;

	// Undefined relative error
	p.rel_err = 0;

	if(p.max_err > tol)
		p.failed = true;

	return p;
}


template<typename T>
prec::estimate_result test_fact(interval k, Real tol, unsigned int n) {

	Real max = 0;
	Real sum = 0;
	Real sum2 = 0;

	T prev = 1;

	for (unsigned int i = k.a; i <= k.b; ++i) {

		T res1 = fact<T>(i);
		T res2 = prev * i;

		real diff;

		// Use a discrete metric to check for equivalence
		if(res1 != res2 ||
		   res1 == 0 || res2 == 0 ||
		   prev > res1)
			diff = 1;
		else
			diff = 0;
		
		if(max < diff)
			max = diff;

		sum += diff;
		sum2 += square(diff);

		prev = res1;
	}

	prec::estimate_result p;
	p.max_err = max;
	p.abs_err = sum;
	p.rms_err = th::sqrt(sum2);
	p.mean_err = sum;

	if(p.abs_err >= 1)
		p.failed = true;

	return p;
}



int main(int argc, char const *argv[]) {

	const real MAX = 1000000;
	const real MIN = -MAX;

	prec::state.outputFolder = "test/";
	prec::state.defaultIterations = 100000;

	prec::setup("core", argc, argv);

		prec::equals("th::square(real)", REAL_LAMBDA(th::square), {
			{1, 1},
			{2, 4},
			{3, 9},
			{0, 0},
			{-1, 1}
		});


		prec::equals("th::cube(real)", REAL_LAMBDA(th::cube), {
			{1, 1},
			{2, 8},
			{3, 27},
			{0, 0},
			{-1, -1}
		});


		prec::estimate(
			"th::sqrt(real)",
			REAL_LAMBDA(th::sqrt),
			REAL_LAMBDA(std::sqrt),
			interval(0, MAX));


		prec::estimate(
			"th::cbrt(real)",
			REAL_LAMBDA(th::cbrt),
			REAL_LAMBDA(std::cbrt),
			interval(-100000, 100000));


		prec::estimate(
			"th::root(real) (2)",
			[](real x) { return th::root(x, 2); },
			REAL_LAMBDA(std::sqrt),
			interval(0, 100000));


		prec::estimate(
			"th::root(real) (3)",
			[](real x) { return th::root(x, 3); },
			REAL_LAMBDA(std::cbrt),
			interval(-100000, 100000));


		prec::estimate(
			"th::root(real) (4)",
			[](real x) { return th::pow(th::root(x, 4), 4); },
			[](real x) { return x; },
			interval(0, 100000));


		prec::estimate(
			"th::isqrt(uint32_t)",
			[](real x) { return th::isqrt<uint32_t>(x); },
			[](real x) { return std::floor(std::sqrt(x)); },
			interval(0, MAX));


		prec::estimate(
			"th::isqrt(uint64_t)",
			[](real x) { return th::isqrt<uint64_t>(x); },
			[](real x) { return std::floor(std::sqrt(x)); },
			interval(0, MAX));


		prec::estimate(
			"th::icbrt(uint32_t)",
			[](real x) { return th::icbrt<uint32_t>(x); },
			[](real x) { return std::floor(std::cbrt(x)); },
			interval(0, MAX));


		prec::estimate(
			"th::icbrt(uint64_t)",
			[](real x) { return th::icbrt<uint64_t>(x); },
			[](real x) { return std::floor(std::cbrt(x)); },
			interval(0, MAX));


		prec::estimate(
			"th::abs(real)",
			REAL_LAMBDA(th::abs),
			REAL_LAMBDA(std::abs),
			interval(MIN, MAX));


		prec::estimate(
			"th::floor(real)",
			REAL_LAMBDA(th::floor),
			REAL_LAMBDA(std::floor),
			interval(MIN, MAX));


		prec::estimate(
			"th::fract(real)",
			REAL_LAMBDA(th::fract),
			[](real x) { return std::abs(std::floor(x) - x); },
			interval(MIN, MAX));


		prec::equals("th::sgn(real)", REAL_LAMBDA(th::sgn), {
				{1, 1},
				{2, 1},
				{-1, -1},
				{-3, -1},
				{0, 0},
				{-1.0 / 3.0, -1}
		});


		prec::estimate(
			"th::ln(real)",
			REAL_LAMBDA(th::ln),
			REAL_LAMBDA(std::log),
			interval(0.00000001, MAX));


		prec::estimate(
			"th::log2(real)",
			REAL_LAMBDA(th::log2),
			REAL_LAMBDA(std::log2),
			interval(0.00000001, MAX));


		prec::estimate(
			"th::log10(real)",
			REAL_LAMBDA(th::log10),
			REAL_LAMBDA(std::log10),
			interval(0.00000001, MAX));


		prec::estimate(
			"th::ilog2(uint32_t)",
			[](real x) { return ilog2<uint32_t>(uint32_t(x)); },
			[](real x) { return uint32_t(std::log2(x)); },
			interval(1, MAX));


		prec::estimate(
			"th::ilog2(uint64_t)",
			[](real x) { return ilog2<uint64_t>(uint64_t(x)); },
			[](real x) { return uint64_t(std::log2(x)); },
			interval(1, MAX));


		prec::estimate(
			"th::pad2(uint32_t)",
			[](real x) { return pad2<uint32_t>(uint32_t(x)); },
			[](real x) { return 1 << (uint32_t) ceil(std::log2(x)); },
			interval(1, MAX));


		prec::estimate(
			"th::pad2(uint64_t)",
			[](real x) { return pad2<uint64_t>(uint64_t(x)); },
			[](real x) { return 1 << (uint64_t) ceil(std::log2(x)); },
			interval(1, MAX));


		prec::estimate(
			"th::exp(real)",
			REAL_LAMBDA(th::exp),
			REAL_LAMBDA(std::exp),
			interval(-100, 10));


		prec::estimate(
			"th::expm1(real)",
			REAL_LAMBDA(th::expm1),
			REAL_LAMBDA(std::expm1),
			interval(-1, 1));


	// test_start("th::pow");

	// 	const unsigned int N = 7;
	// 	const unsigned int MAX_POW = 10;

	// 	for (unsigned int i = 0; i < N; ++i) {
	// 		for (unsigned int j = 0; j < 100; ++j) {
	// 			real x = 20 * i / N;
	// 			real p = MAX_POW * j / 100;
	// 			test_tol(th::pow(x, p), std::pow(x, p), x, TOLERANCE, true);
	// 		}
	// 	}

	// test_end();


		{
			prec::equals("th::powf", th::powf(2, 0.5), th::SQRT2);
			prec::equals("th::powf", th::powf(2, -0.5), 1 / th::SQRT2);
			prec::equals("th::powf", th::powf(2, 2), 4);
			prec::equals("th::powf", th::powf(3, 2), 9);
		}


		prec::estimate(
			"th::sin(real)",
			REAL_LAMBDA(th::sin),
			REAL_LAMBDA(std::sin),
			interval(MIN, MAX));


		prec::estimate(
			"th::cos(real)",
			REAL_LAMBDA(th::cos),
			REAL_LAMBDA(std::cos),
			interval(MIN, MAX));


		prec::estimate(
			"sin^2 + cos^2 = 1",
			[](real x) {
				return square(th::sin(x)) + square(th::cos(x));
			},
			[](real x) {
				return 1;
			},
			interval(MIN, MAX));


		prec::estimate(
			"th::tan(real)",
			REAL_LAMBDA(th::tan),
			REAL_LAMBDA(std::tan),
			interval(-1, 1));


		prec::estimate(
			"th::tan(real)",
			REAL_LAMBDA(th::tan),
			REAL_LAMBDA(std::tan),
			interval(MIN, MAX), 1E-6, false, 1000);

		prec::equals("tan(2) = tan(2 + 100 PI)",
			th::tan(2), th::tan(2 + 100 * PI));


		prec::estimate(
			"th::asin(real)",
			REAL_LAMBDA(th::asin),
			REAL_LAMBDA(std::asin),
			interval(-0.999999, 0.999999),
			0.0001);


		prec::estimate(
			"th::acos(real)",
			REAL_LAMBDA(th::acos),
			REAL_LAMBDA(std::acos),
			interval(-0.999999, 0.999999),
			0.0001);


		prec::estimate(
			"th::atan(real)",
			REAL_LAMBDA(th::atan),
			REAL_LAMBDA(std::atan),
			interval(MIN, MAX),
			0.0001);


		prec::estimate(
			"th::sinh(real)",
			REAL_LAMBDA(th::sinh),
			REAL_LAMBDA(std::sinh),
			interval(-10, 10));


		prec::estimate(
			"th::cosh(real)",
			REAL_LAMBDA(th::cosh),
			REAL_LAMBDA(std::cosh),
			interval(-10, 10));


		prec::estimate(
			"th::tanh(real)",
			REAL_LAMBDA(th::tanh),
			REAL_LAMBDA(std::tanh),
			interval(-10, 10));

		{
			prec::equals("th::binomial_coeff", th::binomial_coeff(1, 1), 1, 0);
			prec::equals("th::binomial_coeff", th::binomial_coeff(2, 0), 1, 0);
			prec::equals("th::binomial_coeff", th::binomial_coeff(2, 1), 2, 0);
			prec::equals("th::binomial_coeff", th::binomial_coeff(3, 2), 3, 0);
			prec::equals("th::binomial_coeff", th::binomial_coeff(3, 1), 3, 0);
			prec::equals("th::binomial_coeff", th::binomial_coeff(6, 3), 20, 0);
			prec::equals("th::binomial_coeff", th::binomial_coeff(10, 3), 120, 0);
			prec::equals("th::binomial_coeff", th::binomial_coeff(16, 7), 11440, 0);
			prec::equals("th::binomial_coeff", th::binomial_coeff(18, 6), 18564, 0);
		}


		prec::equals("th::degrees(real)", REAL_LAMBDA(th::degrees), {
			{th::PI, 180},
			{th::PI/2.0, 90},
			{th::PI/4.0, 45},
			{0, 0},
		});


		prec::equals("th::radians(real)", REAL_LAMBDA(th::radians), {
			{180, th::PI},
			{90, th::PI/2.0},
			{45, th::PI/4.0},
			{0, 0},
		});


		// Square a relatively small number and check that the high bits are zero
		prec::estimate("th::mul_uint128",
			[](real x) {

				uint64_t i = (uint64_t) x;
				uint64_t r1, r2;
				mul_uint128(i, i, r1, r2);

				return r2;
			}, [](real x) { return 0; }, interval(0, 1000));


		prec::estimate("ratio::eval<real>", test_ratio, interval(MIN, MAX));


		prec::estimate("fact<uint32_t>",
			test_fact<uint32_t>, interval(1, 13));


		prec::estimate("fact<uint64_t>",
			test_fact<uint64_t>, interval(1, 20));


		prec::estimate(
			"falling_fact (0)",
			[](real x) {
				return falling_fact(x, 0);
			},
			[](real x) { return 1; },
			interval(0, MAX)
		);


		prec::estimate(
			"falling_fact (1)",
			[](real x) {
				return falling_fact(x, 1);
			},
			[](real x) { return x; },
			interval(0, MAX)
		);


		prec::estimate(
			"falling_fact (2)",
			[](real x) {
				return falling_fact(x, 2);
			},
			[](real x) { return square(x) - x; },
			interval(0, MAX)
		);


		prec::estimate(
			"falling_fact (3)",
			[](real x) {
				return falling_fact(x, 3);
			},
			[](real x) { return cube(x) - 3 * square(x) + 2 * x; },
			interval(0, 100000)
		);


		prec::estimate(
			"rising_fact (0)",
			[](real x) {
				return rising_fact(x, 0);
			},
			[](real x) { return 1; },
			interval(0, MAX)
		);


		prec::estimate(
			"rising_fact (1)",
			[](real x) {
				return rising_fact(x, 1);
			},
			[](real x) { return x; },
			interval(0, MAX)
		);


		prec::estimate(
			"rising_fact (2)",
			[](real x) {
				return rising_fact(x, 2);
			},
			[](real x) { return square(x) + x; },
			interval(0, MAX)
		);


		prec::estimate(
			"rising_fact (3)",
			[](real x) {
				return rising_fact(x, 3);
			},
			[](real x) { return cube(x) + 3 * square(x) + 2 * x; },
			interval(0, 100000)
		);


		// Special functions


		// Gamma function

		// Check translation identity
		prec::estimate(
			"gamma (1)",
			[](real x) {
				return special::gamma(x);
			},
			[](real x) {
				return special::gamma(x + 1) / x;
			},
			interval(0.1, 20),
			10E-8, false,
			prec::state.defaultIterations, prec::fail_on_rel_err
		);

		// Check identity with factorial
		prec::estimate(
			"gamma (2)",
			[](real x) {
				return special::gamma((real) th::floor(x));
			},
			[](real x) {
				return (real) fact<uint64_t>((unsigned int) (th::floor(x) - 1));
			},
			interval(1, 20)
		);


		// Pi function

		// Check identity with factorial
		prec::estimate(
			"pi (2)",
			[](real x) {
				return special::pi((real) th::floor(x));
			},
			[](real x) {
				return (real) fact<uint64_t>((unsigned int) (th::floor(x)));
			},
			interval(1, 20)
		);


	prec::terminate();
}
