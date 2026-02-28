
/// @file test_core.cpp Test cases for real functions and core functionalities

#include "theoretica.h"
#include "chebyshev.h"
#include <ctime>
#include <cmath>

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {


	auto ctx = prec::make_context("core", argc, argv);
	ctx.output->settings.outputFiles = { "test/prec/prec_core.csv" };
	ctx.settings.defaultIterations = 1'000'000;


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

	ctx.estimate(
		"th::sqrt(real)",
		CAST_LAMBDA(th::sqrt, real),
		CAST_LAMBDA(std::sqrt, real),
		Rplus_opt
	);

	ctx.estimate(
		"th::sqrt^2 = th::abs",
		[](real x) { return th::square(th::sqrt(x)); },
		[](real x) { return th::abs(x); },
		Rplus_opt
	);

	ctx.estimate(
		"th::cbrt(real)",
		CAST_LAMBDA(th::cbrt, real),
		CAST_LAMBDA(std::cbrt, real),
		R_opt
	);

	ctx.estimate(
		"th::cbrt^3(x) = x",
		[](real x) { return th::cube(th::cbrt(x)); },
		[](real x) { return x; },
		R_opt
	);

	ctx.estimate(
		"th::root(real) (2)",
		[](real x) { return th::root(x, 2); },
		CAST_LAMBDA(std::sqrt, real),
		Rplus_opt
	);


	ctx.estimate(
		"th::root(real) (3)",
		[](real x) { return th::root(x, 3); },
		CAST_LAMBDA(std::cbrt, real),
		R_opt
	);


	ctx.estimate(
		"th::root(real) (4)",
		[](real x) { return th::pow(th::root(x, 4), 4); },
		[](real x) { return x; },
		Rplus_opt
	);


	ctx.estimate(
		"th::isqrt(uint32_t)",
		th::isqrt<uint32_t>,
		[](real x) { return std::floor(std::sqrt(x)); },
		Rplus_opt
	);


	ctx.estimate(
		"th::isqrt(uint64_t)",
		th::isqrt<uint64_t>,
		[](real x) { return std::floor(std::sqrt(x)); },
		Rplus_opt
	);


	ctx.estimate(
		"th::icbrt(uint32_t)",
		th::icbrt<uint32_t>,
		[](real x) { return std::floor(std::cbrt(x)); },
		Rplus_opt
	);


	ctx.estimate(
		"th::icbrt(uint64_t)",
		th::icbrt<uint64_t>,
		[](real x) { return std::floor(std::cbrt(x)); },
		Rplus_opt
	);


	ctx.estimate(
		"th::abs(real)",
		CAST_LAMBDA(th::abs, real),
		CAST_LAMBDA(std::abs, real),
		R_opt
	);


	ctx.estimate(
		"th::floor(real)",
		CAST_LAMBDA(th::floor, real),
		CAST_LAMBDA(std::floor, real),
		R_opt
	);


	ctx.estimate(
		"th::fract(real)",
		CAST_LAMBDA(th::fract, real),
		[](real x) { return std::abs(std::floor(x) - x); },
		R_opt
	);


	ctx.estimate(
		"th::sgn (1)",
		CAST_LAMBDA(th::sgn, real),
		[](real x) { return 1; },
		prec::estimate_options<real, real>(
			prec::interval(0.1, 1E+06),
			prec::estimator::quadrature1D()
		)
	);


	ctx.estimate(
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


	ctx.estimate(
		"th::ln(real)",
		CAST_LAMBDA(th::ln, real),
		CAST_LAMBDA(std::log, real),
		log_opt
	);


	ctx.estimate(
		"th::log2(real)",
		CAST_LAMBDA(th::log2, real),
		CAST_LAMBDA(std::log2, real),
		log_opt
	);
	
	ctx.estimate(
		"th::log10(real)",
		CAST_LAMBDA(th::log10, real),
		CAST_LAMBDA(std::log10, real),
		log_opt
	);


	ctx.estimate(
		"th::ilog2(uint32_t)",
		[](real x) { return ilog2<uint32_t>(uint32_t(x)); },
		[](real x) { return uint32_t(std::log2(x)); },
		prec::interval(1, 1E+06)
	);


	ctx.estimate(
		"th::ilog2(uint64_t)",
		[](real x) { return ilog2<uint64_t>(uint64_t(x)); },
		[](real x) { return uint64_t(std::log2(x)); },
		prec::interval(1, 1E+06)
	);


	ctx.estimate(
		"th::pad2(uint32_t)",
		[](real x) { return pad2<uint32_t>(uint32_t(x)); },
		[](real x) { return 1 << (uint32_t) ceil(std::log2(x)); },
		prec::estimate_options<real, real>(
			prec::interval(1, 1E+06),
			prec::estimator::quadrature1D()
		)
	);


	ctx.estimate(
		"th::pad2(uint64_t)",
		[](real x) { return pad2<uint64_t>(uint64_t(x)); },
		[](real x) { return 1 << (uint64_t) ceil(std::log2(x)); },
		prec::estimate_options<real, real>(
			prec::interval(1, 1E+06),
			prec::estimator::quadrature1D()
		)
	);


	ctx.estimate(
		"th::exp(real)",
		CAST_LAMBDA(th::exp, real),
		CAST_LAMBDA(std::exp, real),
		exp_opt
	);


	ctx.estimate(
		"th::expm1(real)",
		CAST_LAMBDA(th::expm1, real),
		CAST_LAMBDA(std::expm1, real),
		prec::estimate_options<real, real>(
			prec::interval(-1, 1),
			prec::estimator::quadrature1D()
		)
	);

	{
		ctx.equals("th::pow", th::pow(1, 1E+06), 1);
		ctx.equals("th::pow", th::pow(1, -1E+06), 1);
		ctx.equals("th::pow", th::pow(2, 10), (1 << 10));
		ctx.equals("th::pow", th::pow(10, 6), 1E+6);
		ctx.equals("th::pow", th::pow(E, 10) * th::pow(E, -10), 1);
		ctx.equals("th::pow", th::pow(1E-08, 10) * th::pow(1E-08, -10), 1);
	}


	{
		ctx.equals("th::ipow", th::ipow(1, 1E+06), 1);
		ctx.equals("th::ipow", th::ipow(2, 10), (1 << 10));
		ctx.equals("th::ipow", th::ipow(10, 6), 1E+6);
	}

	{
		ctx.equals("th::powf", th::powf(2, 0.5), th::SQRT2);
		ctx.equals("th::powf", th::powf(2, -0.5), 1 / th::SQRT2);
		ctx.equals("th::powf", th::powf(2, 2), 4);
		ctx.equals("th::powf", th::powf(3, 2), 9);
	}


	ctx.estimate(
		"th::sin(real)",
		CAST_LAMBDA(th::sin, real),
		CAST_LAMBDA(std::sin, real),
		R_opt
	);


	ctx.estimate(
		"th::cos(real)",
		CAST_LAMBDA(th::cos, real),
		CAST_LAMBDA(std::cos, real),
		R_opt
	);


	ctx.estimate(
		"sin^2 + cos^2 = 1",
		[](real x) { return square(th::sin(x)) + square(th::cos(x)); },
		[](real x) { return 1; },
		R_opt
	);


	ctx.estimate(
		"th::tan(real)",
		CAST_LAMBDA(th::tan, real),
		CAST_LAMBDA(std::tan, real),
		prec::estimate_options<real, real>(
			prec::interval(-1, 1),
			prec::estimator::quadrature1D()
		)
	);

	ctx.equals(
		"tan(2)=tan(2+100*PI)",
		th::tan(2.0),
		th::tan(2.0 + 100 * PI)
	);

	auto asin_opt = prec::estimate_options<real, real>(
		prec::interval(-0.999999, +0.999999),
		prec::estimator::quadrature1D()
	);

	ctx.estimate(
		"th::asin(real)",
		CAST_LAMBDA(th::asin, real),
		CAST_LAMBDA(std::asin, real),
		asin_opt	
	);


	ctx.estimate(
		"th::acos(real)",
		CAST_LAMBDA(th::acos, real),
		CAST_LAMBDA(std::acos, real),
		asin_opt
	);


	ctx.estimate(
		"th::atan(real)",
		CAST_LAMBDA(th::atan, real),
		CAST_LAMBDA(std::atan, real),
		R_opt
	);


	ctx.estimate(
		"th::sinh(real)",
		CAST_LAMBDA(th::sinh, real),
		CAST_LAMBDA(std::sinh, real),
		exp_opt
	);


	ctx.estimate(
		"th::cosh(real)",
		CAST_LAMBDA(th::cosh, real),
		CAST_LAMBDA(std::cosh, real),
		exp_opt
	);


	ctx.estimate(
		"th::tanh(real)",
		CAST_LAMBDA(th::tanh, real),
		CAST_LAMBDA(std::tanh, real),
		exp_opt
	);

	{
		ctx.equals("th::binomial_coeff", th::binomial_coeff(1, 1), 1, 0);
		ctx.equals("th::binomial_coeff", th::binomial_coeff(2, 0), 1, 0);
		ctx.equals("th::binomial_coeff", th::binomial_coeff(2, 1), 2, 0);
		ctx.equals("th::binomial_coeff", th::binomial_coeff(3, 2), 3, 0);
		ctx.equals("th::binomial_coeff", th::binomial_coeff(3, 1), 3, 0);
		ctx.equals("th::binomial_coeff", th::binomial_coeff(6, 3), 20, 0);
		ctx.equals("th::binomial_coeff", th::binomial_coeff(10, 3), 120, 0);
		ctx.equals("th::binomial_coeff", th::binomial_coeff(16, 7), 11440, 0);
		ctx.equals("th::binomial_coeff", th::binomial_coeff(18, 6), 18564, 0);
	}


// 	ctx.equals("th::degrees(real)", CAST_LAMBDA(th::degrees), {
// 		{th::PI, 180},
// 		{th::PI/2.0, 90},
// 		{th::PI/4.0, 45},
// 		{0, 0},
// 	});


// 	ctx.equals("th::radians(real)", CAST_LAMBDA(th::radians), {
// 		{180, th::PI},
// 		{90, th::PI/2.0},
// 		{45, th::PI/4.0},
// 		{0, 0},
// 	});


	// Square a relatively small number and check that the high bits are zero
	ctx.estimate("th::mul_uint128",
		[](real x) {

			uint64_t i = (uint64_t) x;
			uint64_t r1, r2;
			mul_uint128(i, i, r1, r2);
			
			return r2;

		}, [](real x) { return 0; }, prec::interval(0, 1000)
	);


// 	ctx.estimate("ratio::eval<real>", test_ratio, R_opt);


	auto fact_opt = prec::estimate_options<uint32_t, uint64_t>(
		prec::interval(1, 20),
		prec::estimator::discrete1D<uint32_t, uint64_t>()
	);

	ctx.estimate("fact<uint32_t>",
		th::fact<uint32_t>,
		[](uint32_t x) { return th::fact<uint32_t>(x - 1) * x; },
		fact_opt
	);


	ctx.estimate("fact<uint64_t>",
		th::fact<uint64_t>,
		[](uint32_t x) { return th::fact<uint64_t>(x - 1) * x; },
		fact_opt
	);


	ctx.estimate(
		"falling_fact(x, 0)",
		[](uint32_t x) { return falling_fact(x, 0); },
		[](uint32_t x) { return 1; },
		fact_opt
	);


	ctx.estimate(
		"falling_fact(x, 1)",
		[](uint32_t x) { return falling_fact(x, 1); },
		[](uint32_t x) { return x; },
		fact_opt
	);


	ctx.estimate(
		"falling_fact(x, 2)",
		[](uint32_t x) { return falling_fact(x, 2); },
		[](uint32_t x) { return square(x) - x; },
		fact_opt
	);


	ctx.estimate(
		"falling_fact(x, 3)",
		[](uint32_t x) { return falling_fact(x, 3); },
		[](uint32_t x) { return cube(x) - 3 * square(x) + 2 * x; },
		fact_opt
	);


	ctx.estimate(
		"rising_fact(x, 0)",
		[](uint32_t x) { return rising_fact(x, 0); },
		[](uint32_t x) { return 1; },
		fact_opt
	);


	ctx.estimate(
		"rising_fact(x, 1)",
		[](uint32_t x) { return rising_fact(x, 1); },
		[](uint32_t x) { return x; },
		fact_opt
	);


	ctx.estimate(
		"rising_fact(x, 2)",
		[](uint32_t x) { return rising_fact(x, 2); },
		[](uint32_t x) { return square(x) + x; },
		fact_opt
	);


	ctx.estimate(
		"rising_fact(x, 3)",
		[](uint32_t x) { return rising_fact(x, 3); },
		[](uint32_t x) { return cube(x) + 3 * square(x) + 2 * x; },
		fact_opt
	);


	// Test special.h

	// Special functions
	auto special_opt = prec::estimate_options<real, real>(
		prec::interval(0.1, 20),
		prec::estimator::quadrature1D()
	);

	special_opt.fail = prec::fail::fail_on_rel_err();

	{
		ctx.equals("special::gamma(uint) gamma(1)", special::gamma(1u), real(1), 0);
		ctx.equals("special::gamma(uint) gamma(2)", special::gamma(2u), real(1), 0);
		ctx.equals("special::gamma(uint) gamma(3)", special::gamma(3u), real(2), 0);
		ctx.equals("special::gamma(uint) gamma(6)", special::gamma(6u), real(120), 0);
		ctx.equals("special::gamma(uint) gamma(0) is NaN", th::is_nan(special::gamma(0u)), true, 0);
	}

	{	
		ctx.equals("special::half_gamma(2)", special::half_gamma(2u), real(1), 1E-12);
		ctx.equals("special::half_gamma(4)", special::half_gamma(4u), real(1), 1E-12);
		ctx.equals("special::half_gamma(6)", special::half_gamma(6u), real(2), 1E-12);

		ctx.equals("special::half_gamma(1)", special::half_gamma(1u), SQRTPI, 1E-8);
		ctx.equals("special::half_gamma(3)", special::half_gamma(3u), SQRTPI / 2.0, 1E-8);
		ctx.equals("special::half_gamma(5)", special::half_gamma(5u), 3.0 * SQRTPI / 4.0, 1E-8);

		ctx.equals("special::half_gamma(0) is NaN", th::is_nan(special::half_gamma(0u)), true, 0);
	}

	ctx.estimate(
		"special::lngamma(real) vs std::lgamma",
		CAST_LAMBDA(special::lngamma, real),
		CAST_LAMBDA(std::lgamma, real),
		special_opt
	);

	{
		ctx.equals(
			"special::lngamma(0.5)",
			special::lngamma(0.5),
			th::ln(SQRTPI),
			1E-8
		);

		ctx.equals(
			"special::lngamma(-1.5)",
			special::lngamma(-1.5),
			real(std::lgamma(-1.5)),
			1E-8
		);

		ctx.equals(
			"special::lngamma(-0.5) is NaN",
			th::is_nan(special::lngamma(-0.5)),
			true,
			0
		);
	}

	ctx.estimate(
		"special::gamma(real) vs std::tgamma",
		CAST_LAMBDA(special::gamma, real),
		CAST_LAMBDA(std::tgamma, real),
		special_opt
	);

	{
		ctx.equals("special::gamma(0.5)", special::gamma(0.5), SQRTPI, 1E-8);
		ctx.equals("special::gamma(1.0)", special::gamma(1.0), real(1), 1E-12);
		ctx.equals("special::gamma(1.5)", special::gamma(1.5), SQRTPI / 2.0, 1E-8);
		ctx.equals("special::gamma(2.5)", special::gamma(2.5), 3.0 * SQRTPI / 4.0, 1E-8);
		ctx.equals("special::gamma(-0.5)", special::gamma(-0.5), -2.0 * SQRTPI, 1E-8);

		ctx.equals("special::gamma(0.0) is inf", th::is_inf(special::gamma(0.0)), true, 0);
		ctx.equals("special::gamma(-1.0) is inf", th::is_inf(special::gamma(-1.0)), true, 0);
		ctx.equals("special::gamma(-2.0) is inf", th::is_inf(special::gamma(-2.0)), true, 0);
	}

	{
		// Integer identity: Pi(n) = n!
		ctx.equals("special::pi(0)", special::pi(0.0), real(1), 1E-12);
		ctx.equals("special::pi(1)", special::pi(1.0), real(1), 1E-12);
		ctx.equals("special::pi(2)", special::pi(2.0), real(2), 1E-12);
		ctx.equals("special::pi(5)", special::pi(5.0), real(120), 1E-8);

		// Relation Pi(x) = Gamma(x + 1)
		ctx.equals(
			"special::pi(0.5) = gamma(1.5)",
			special::pi(0.5),
			special::gamma(1.5),
			1E-8
		);
	}

	{
		ctx.equals("special::beta(1,1)", special::beta(1.0, 1.0), real(1), 1E-12);
		ctx.equals("special::beta(1,2)", special::beta(1.0, 2.0), real(0.5), 1E-10);
		ctx.equals("special::beta(2,3)", special::beta(2.0, 3.0), real(1.0 / 12.0), 1E-10);
		ctx.equals("special::beta(0.5,0.5)", special::beta(0.5, 0.5), PI, 1E-7);

		ctx.equals("special::beta symmetry (0.5, 1.5)", special::beta(0.5, 1.5), special::beta(1.5, 0.5), 1E-8);
		ctx.equals("special::beta symmetry (1, 3)", special::beta(1.0, 3.0), special::beta(3.0, 1.0), 1E-8);
		ctx.equals("special::beta symmetry (2.5, 4)", special::beta(2.5, 4.0), special::beta(4.0, 2.5), 1E-8);

	}


	// Test bit_op.h
	

	{
		uint64_t a = 0xFFFFFFFFFFFFFFFF;
		uint64_t b = 0x2;
		uint64_t c_low, c_high;

		uint64_t expected_low = 0xFFFFFFFFFFFFFFFE;
		uint64_t expected_high = 0x1;

		th::mul_uint128(a, b, c_low, c_high);

		ctx.equals("th::mul_uint128 (c_low)", c_low, expected_low);
		ctx.equals("th::mul_uint128 (c_high)", c_high, expected_high);
	}

	{
		uint64_t a = 0;
		uint64_t b = 0;
		
		uint64_t result = th::mix_mum(a, b);
		
		ctx.equals("th::mix_mum == 0", result, 0);
	}

	{
		uint64_t a = 0x12345678ABCDEF00;
		uint64_t b = 0x0FEDCBA987654321;
		
		uint64_t result = th::mix_mum(a, b);
		
		ctx.equals("th::mix_mum != 0", result != 0, true, 0);
	}

	{
		uint64_t x = 0x12345678ABCDEF00;
		unsigned int i = 8;
		uint64_t rotated = th::bit_rotate(x, i);

		uint64_t expected_rotated = 0x345678ABCDEF0012;

		ctx.equals("th::bit_rotate (64-bit)", rotated, expected_rotated);
	}

	{
		uint32_t x = 0xABCDEF00;
		unsigned int i = 4;
		uint32_t rotated = th::bit_rotate(x, i);

		uint32_t expected_rotated = 0xBCDEF00A;

		ctx.equals("th::bit_rotate (32-bit)", rotated, expected_rotated);
	}

	{
		std::vector<uint8_t> v = {};
		
		th::swap_bit_reverse(v, 0);

		std::vector<uint8_t> expected = {};

		ctx.equals("th::swap_bit_reverse", v == expected, true);
	}

	{
		std::vector<uint8_t> v = {1};
		
		th::swap_bit_reverse(v, 0);

		std::vector<uint8_t> expected = {1};

		ctx.equals("th::swap_bit_reverse", v == expected, true);
	}

	{
		std::vector<uint8_t> v = {1, 2, 3, 4};
		
		th::swap_bit_reverse(v, 2);

		std::vector<uint8_t> expected = {1, 3, 2, 4};

		ctx.equals("th::swap_bit_reverse", v == expected, true);
	}

	{
		vec<uint8_t> v = {1, 2, 3, 4, 5, 6};
		
		th::swap_bit_reverse(v, 2);

		vec<uint8_t> expected = {1, 3, 2, 4, 5, 6};

		ctx.equals("th::swap_bit_reverse", v == expected, true);
	}


	// error.h

	ctx.equals("th::is_nan", th::is_nan(1.0), false);
	ctx.equals("th::is_nan", th::is_nan(th::nan()), true);

	ctx.equals("th::is_inf", th::is_inf(1.0), false);
	ctx.equals("th::is_inf", th::is_inf(th::inf()), true);


	// reprod.h
	auto env = reprod::get_env();
	ctx.equals("get_env().os", env.os != "", true);
	ctx.equals("get_env().arch", env.arch != "", true);
	ctx.equals("get_env().compiler", env.compiler != "", true);
	ctx.equals("get_env().compiler_version", env.compiler_version != "", true);
	ctx.equals("get_env().build_date", env.build_date != "", true);
	ctx.equals("get_env().cpp_standard", env.cpp_standard != "", true);
}
