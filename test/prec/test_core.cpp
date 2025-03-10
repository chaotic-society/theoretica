
/// @file test_core.cpp Test cases for real functions and core functionalities

#include "theoretica.h"
#include "chebyshev.h"
#include <ctime>
#include <cmath>

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {

	prec::setup("core", argc, argv);

		output::settings.outputFiles = { "test/prec/prec_core.csv" };

		prec::settings.defaultIterations = 1'000'000;

		// Test real_analysis.h


		// Estimate options for real endofunctions.
		auto R_opt = prec::estimate_options<real, real>(
				prec::interval(-1E+06, 1E+06),
				prec::estimator::quadrature1D()
		);

		// Estimate options for functions defined
		// over the positive real numbers.
		auto Rplus_opt = prec::estimate_options<real, real>(
				prec::interval(0, 1E+06),
				prec::estimator::quadrature1D()
		);

		auto exp_opt = prec::estimate_options<real, real>(
			prec::interval(-10, 10),
			prec::estimator::quadrature1D()
		);

		prec::estimate(
			"th::sqrt(real)",
			CAST_LAMBDA(th::sqrt, real),
			CAST_LAMBDA(std::sqrt, real),
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
			CAST_LAMBDA(th::cbrt, real),
			CAST_LAMBDA(std::cbrt, real),
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
			CAST_LAMBDA(std::sqrt, real),
			Rplus_opt
		);


		prec::estimate(
			"th::root(real) (3)",
			[](real x) { return th::root(x, 3); },
			CAST_LAMBDA(std::cbrt, real),
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
			CAST_LAMBDA(th::abs, real),
			CAST_LAMBDA(std::abs, real),
			R_opt
		);


		prec::estimate(
			"th::floor(real)",
			CAST_LAMBDA(th::floor, real),
			CAST_LAMBDA(std::floor, real),
			R_opt
		);


		prec::estimate(
			"th::fract(real)",
			CAST_LAMBDA(th::fract, real),
			[](real x) { return std::abs(std::floor(x) - x); },
			R_opt
		);


		prec::estimate(
			"th::sgn (1)",
			CAST_LAMBDA(th::sgn, real),
			[](real x) { return 1; },
			prec::estimate_options<real, real>(
				prec::interval(0.1, 1E+06),
				prec::estimator::quadrature1D()
			)
		);


		prec::estimate(
			"th::sgn (2)",
			CAST_LAMBDA(th::sgn, real),
			[](real x) { return -1; },
			prec::estimate_options<real, real>(
				prec::interval(-1E+06, -0.1),
				prec::estimator::quadrature1D()
			)
		);


		auto log_opt = prec::estimate_options<real, real>(
			prec::interval(1E-08, 1E+06),
			prec::estimator::quadrature1D()
		);


		prec::estimate(
			"th::ln(real)",
			CAST_LAMBDA(th::ln, real),
			CAST_LAMBDA(std::log, real),
			log_opt
		);


		prec::estimate(
			"th::log2(real)",
			CAST_LAMBDA(th::log2, real),
			CAST_LAMBDA(std::log2, real),
			log_opt
		);


		prec::estimate(
			"th::log10(real)",
			CAST_LAMBDA(th::log10, real),
			CAST_LAMBDA(std::log10, real),
			log_opt
		);


		prec::estimate(
			"th::ilog2(uint32_t)",
			[](real x) { return ilog2<uint32_t>(uint32_t(x)); },
			[](real x) { return uint32_t(std::log2(x)); },
			Rplus_opt
		);


		prec::estimate(
			"th::ilog2(uint64_t)",
			[](real x) { return ilog2<uint64_t>(uint64_t(x)); },
			[](real x) { return uint64_t(std::log2(x)); },
			Rplus_opt
		);


		prec::estimate(
			"th::pad2(uint32_t)",
			[](real x) { return pad2<uint32_t>(uint32_t(x)); },
			[](real x) { return 1 << (uint32_t) ceil(std::log2(x)); },
			prec::estimate_options<real, real>(
				prec::interval(1, 1E+06),
				prec::estimator::quadrature1D()
			)
		);


		prec::estimate(
			"th::pad2(uint64_t)",
			[](real x) { return pad2<uint64_t>(uint64_t(x)); },
			[](real x) { return 1 << (uint64_t) ceil(std::log2(x)); },
			prec::estimate_options<real, real>(
				prec::interval(1, 1E+06),
				prec::estimator::quadrature1D()
			)
		);


		prec::estimate(
			"th::exp(real)",
			CAST_LAMBDA(th::exp, real),
			CAST_LAMBDA(std::exp, real),
			exp_opt
		);


		prec::estimate(
			"th::expm1(real)",
			CAST_LAMBDA(th::expm1, real),
			CAST_LAMBDA(std::expm1, real),
			prec::estimate_options<real, real>(
				prec::interval(-1, 1),
				prec::estimator::quadrature1D()
			)
		);

		{
			prec::equals("th::pow", th::pow(1, 1E+06), 1);
			prec::equals("th::pow", th::pow(1, -1E+06), 1);
			prec::equals("th::pow", th::pow(2, 10), (1 << 10));
			prec::equals("th::pow", th::pow(10, 6), 1E+6);
			prec::equals("th::pow", th::pow(E, 10) * th::pow(E, -10), 1);
			prec::equals("th::pow", th::pow(1E-08, 10) * th::pow(1E-08, -10), 1);
		}


		{
			prec::equals("th::ipow", th::ipow(1, 1E+06), 1);
			prec::equals("th::ipow", th::ipow(2, 10), (1 << 10));
			prec::equals("th::ipow", th::ipow(10, 6), 1E+6);
		}

		{
			prec::equals("th::powf", th::powf(2, 0.5), th::SQRT2);
			prec::equals("th::powf", th::powf(2, -0.5), 1 / th::SQRT2);
			prec::equals("th::powf", th::powf(2, 2), 4);
			prec::equals("th::powf", th::powf(3, 2), 9);
		}


		prec::estimate(
			"th::sin(real)",
			CAST_LAMBDA(th::sin, real),
			CAST_LAMBDA(std::sin, real),
			R_opt
		);


		prec::estimate(
			"th::cos(real)",
			CAST_LAMBDA(th::cos, real),
			CAST_LAMBDA(std::cos, real),
			R_opt
		);


		prec::estimate(
			"sin^2 + cos^2 = 1",
			[](real x) { return square(th::sin(x)) + square(th::cos(x)); },
			[](real x) { return 1; },
			R_opt
		);


		prec::estimate(
			"th::tan(real)",
			CAST_LAMBDA(th::tan, real),
			CAST_LAMBDA(std::tan, real),
			prec::estimate_options<real, real>(
				prec::interval(-1, 1),
				prec::estimator::quadrature1D()
			)
		);

		prec::equals(
			"tan(2)=tan(2+100*PI)",
			th::tan(2.0),
			th::tan(2.0 + 100 * PI)
		);

		auto asin_opt = prec::estimate_options<real, real>(
			prec::interval(-0.999999, +0.999999),
			prec::estimator::quadrature1D()
		);

		prec::estimate(
			"th::asin(real)",
			CAST_LAMBDA(th::asin, real),
			CAST_LAMBDA(std::asin, real),
			asin_opt	
		);


		prec::estimate(
			"th::acos(real)",
			CAST_LAMBDA(th::acos, real),
			CAST_LAMBDA(std::acos, real),
			asin_opt
		);


		prec::estimate(
			"th::atan(real)",
			CAST_LAMBDA(th::atan, real),
			CAST_LAMBDA(std::atan, real),
			R_opt
		);


		prec::estimate(
			"th::sinh(real)",
			CAST_LAMBDA(th::sinh, real),
			CAST_LAMBDA(std::sinh, real),
			exp_opt
		);


		prec::estimate(
			"th::cosh(real)",
			CAST_LAMBDA(th::cosh, real),
			CAST_LAMBDA(std::cosh, real),
			exp_opt
		);


		prec::estimate(
			"th::tanh(real)",
			CAST_LAMBDA(th::tanh, real),
			CAST_LAMBDA(std::tanh, real),
			exp_opt
		);

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


	// 	prec::equals("th::degrees(real)", CAST_LAMBDA(th::degrees), {
	// 		{th::PI, 180},
	// 		{th::PI/2.0, 90},
	// 		{th::PI/4.0, 45},
	// 		{0, 0},
	// 	});


	// 	prec::equals("th::radians(real)", CAST_LAMBDA(th::radians), {
	// 		{180, th::PI},
	// 		{90, th::PI/2.0},
	// 		{45, th::PI/4.0},
	// 		{0, 0},
	// 	});


		// Square a relatively small number and check that the high bits are zero
		// prec::estimate("th::mul_uint128",
		// 	[](real x) {

		// 		uint64_t i = (uint64_t) x;
		// 		uint64_t r1, r2;
		// 		mul_uint128(i, i, r1, r2);

		// 		return r2;
		// 	}, [](real x) { return 0; }, interval(0, 1000));


	// 	prec::estimate("ratio::eval<real>", test_ratio, R_opt);


		auto fact_opt = prec::estimate_options<uint32_t, uint64_t>(
			prec::interval(1, 20),
			prec::estimator::discrete1D<uint32_t, uint64_t>()
		);

		prec::estimate("fact<uint32_t>",
			th::fact<uint32_t>,
			[](uint32_t x) { return th::fact<uint32_t>(x - 1) * x; },
			fact_opt
		);


		prec::estimate("fact<uint64_t>",
			th::fact<uint64_t>,
			[](uint32_t x) { return th::fact<uint64_t>(x - 1) * x; },
			fact_opt
		);


		prec::estimate(
			"falling_fact(x, 0)",
			[](uint32_t x) { return falling_fact(x, 0); },
			[](uint32_t x) { return 1; },
			fact_opt
		);


		prec::estimate(
			"falling_fact(x, 1)",
			[](uint32_t x) { return falling_fact(x, 1); },
			[](uint32_t x) { return x; },
			fact_opt
		);


		prec::estimate(
			"falling_fact(x, 2)",
			[](uint32_t x) { return falling_fact(x, 2); },
			[](uint32_t x) { return square(x) - x; },
			fact_opt
		);


		prec::estimate(
			"falling_fact(x, 3)",
			[](uint32_t x) { return falling_fact(x, 3); },
			[](uint32_t x) { return cube(x) - 3 * square(x) + 2 * x; },
			fact_opt
		);


		prec::estimate(
			"rising_fact(x, 0)",
			[](uint32_t x) { return rising_fact(x, 0); },
			[](uint32_t x) { return 1; },
			fact_opt
		);


		prec::estimate(
			"rising_fact(x, 1)",
			[](uint32_t x) { return rising_fact(x, 1); },
			[](uint32_t x) { return x; },
			fact_opt
		);


		prec::estimate(
			"rising_fact(x, 2)",
			[](uint32_t x) { return rising_fact(x, 2); },
			[](uint32_t x) { return square(x) + x; },
			fact_opt
		);


		prec::estimate(
			"rising_fact(x, 3)",
			[](uint32_t x) { return rising_fact(x, 3); },
			[](uint32_t x) { return cube(x) + 3 * square(x) + 2 * x; },
			fact_opt
		);


		// Test special.h


		// Special functions
		auto special_opt = prec::estimate_options<real, real>(
			prec::interval(1, 20),
			prec::estimator::quadrature1D()
		);

		special_opt.fail = prec::fail::fail_on_rel_err();


		// Check translation identity
		prec::estimate(
			"special::gamma (1)",
			CAST_LAMBDA(special::gamma, real),
			[](real x) { return special::gamma(x + 1) / x; },
			special_opt
		);


		// Check identity with factorial
		prec::estimate(
			"special::gamma (2)",
			[](real x) { return special::gamma((real) th::floor(x)); },
			[](real x) { return (real) fact<uint64_t>((unsigned int) (th::floor(x) - 1)); },
			special_opt
		);


		// Check identity with factorial
		prec::estimate(
			"special::pi",
			[](real x) {
				return special::pi((real) th::floor(x));
			},
			[](real x) {
				return (real) fact<uint64_t>((unsigned int) (th::floor(x)));
			},
			special_opt
		);

		prec::estimate(
		    "special::half_gamma(uint32_t)",
		    CAST_LAMBDA(special::half_gamma, uint32_t),
		    [](uint32_t k) { 
		        return (k % 2 == 0)
		            ? fact<uint64_t>(k / 2 - 1)
		            : double_fact<uint64_t>(k - 2) * SQRTPI / (1 << ((k - 1) / 2));
		    },
		    special_opt
		);

		prec::estimate(
		    "special::lngamma(real)",
		    CAST_LAMBDA(special::lngamma, real),
		    [](real x) {
		        const real c[7] = {1.000000000178, 76.180091729400, -86.505320327112, 24.014098222230, -1.231739516140, 0.001208580030, -0.000005363820};
		        real A5 = c[0];
		        for (int i = 1; i < 7; ++i)
		            A5 += c[i] / (x + i - 1);
		        return (x - 0.5) * (std::log(x + 4.5) - 1) - 5 + std::log(SQRTPI * SQRT2 * A5);
		    },
		    special_opt
		);

		prec::equals(
		    "special::beta(real)",
		    special::beta(2.0, 3.0),
		    std::exp(special::lngamma(2.0) + special::lngamma(3.0) - special::lngamma(5.0)),
		    1e-8
		);


		// Test bit_op.h
	

	{
		uint64_t a = 0xFFFFFFFFFFFFFFFF;
		uint64_t b = 0x2;
		uint64_t c_low, c_high;

		uint64_t expected_low = 0xFFFFFFFFFFFFFFFE;
		uint64_t expected_high = 0x1;

		th::mul_uint128(a, b, c_low, c_high);

		prec::equals("th::mul_uint128 (c_low)", c_low, expected_low);
		prec::equals("th::mul_uint128 (c_high)", c_high, expected_high);
	}

	{
		uint64_t a = 0;
		uint64_t b = 0;
		
		uint64_t result = th::mix_mum(a, b);
		
		prec::equals("th::mix_mum == 0", result, 0);
	}

	{
		uint64_t a = 0x12345678ABCDEF00;
		uint64_t b = 0x0FEDCBA987654321;
		
		uint64_t result = th::mix_mum(a, b);
		
		prec::equals("th::mix_mum != 0", result != 0, true, 0);
	}

	{
		uint64_t x = 0x12345678ABCDEF00;
		unsigned int i = 8;
		uint64_t rotated = th::bit_rotate(x, i);

		uint64_t expected_rotated = 0x345678ABCDEF0012;

		prec::equals("th::bit_rotate (64-bit)", rotated, expected_rotated);
	}

	{
		uint32_t x = 0xABCDEF00;
		unsigned int i = 4;
		uint32_t rotated = th::bit_rotate(x, i);

		uint32_t expected_rotated = 0xBCDEF00A;

		prec::equals("th::bit_rotate (32-bit)", rotated, expected_rotated);
	}

	{
		std::vector<uint8_t> v = {};
		
		th::swap_bit_reverse(v, 0);

		std::vector<uint8_t> expected = {};

		prec::equals("th::swap_bit_reverse", v == expected, true);
	}

	{
		std::vector<uint8_t> v = {1};
		
		th::swap_bit_reverse(v, 0);

		std::vector<uint8_t> expected = {1};

		prec::equals("th::swap_bit_reverse", v == expected, true);
	}

	{
		std::vector<uint8_t> v = {1, 2, 3, 4};
		
		th::swap_bit_reverse(v, 2);

		std::vector<uint8_t> expected = {1, 3, 2, 4};

		prec::equals("th::swap_bit_reverse", v == expected, true);
	}

	{
		vec<uint8_t> v = {1, 2, 3, 4, 5, 6};
		
		th::swap_bit_reverse(v, 2);

		vec<uint8_t> expected = {1, 3, 2, 4, 5, 6};

		prec::equals("th::swap_bit_reverse", v == expected, true);
	}

	prec::terminate();
}