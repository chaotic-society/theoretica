
///
/// @file dual2_functions.h Functions defined on second order
/// dual numbers for automatic differentiation.
///


#ifndef THEORETICA_DUAL2_FUNCTIONS_H
#define THEORETICA_DUAL2_FUNCTIONS_H

#include "./dual2.h"
#include "../core/real_analysis.h"

#include <iostream>


namespace theoretica {


	// Second derivative of the composite function:
	// (f(g))'' = f''(g) * g'' + f'(g) * g''


	/// Return the square of a second order dual number
	inline dual2 square(dual2 x) {
		return x * x;
	}


	/// Return the cube of a second order dual number
	inline dual2 cube(dual2 x) {
		return x * x * x;
	}


	/// Compute the n-th power of a second order dual number
	inline dual2 pow(dual2 x, int n) {

		real pow_n_2_x = pow(x.Re(), n - 2);

		real f = pow_n_2_x * square(x.Re());
		real df = pow_n_2_x * x.Re() * n * x.Dual1();
		real d2f = (pow_n_2_x * n * (n - 1)) * square(x.Dual1()) + df * x.Dual2();

		return dual2(f, df, d2f);
	}


	/// Compute the square root of a second order dual number
	inline dual2 sqrt(dual2 x) {

		real sqrt_x = sqrt(x.Re());

		if(sqrt_x == 0) {
			TH_MATH_ERROR("sqrt(dual2)", sqrt_x, DIV_BY_ZERO);
			return dual2(nan(), nan());
		}

		return dual2(sqrt_x,
			0.5 / sqrt_x * x.Dual1(),
			-0.25 / x.Re() / sqrt_x * square(x.Dual1()) + 0.5 / sqrt_x * x.Dual2());
	}


	/// Compute the sine of a second order dual number
	inline dual2 sin(dual2 x) {

		real sin_x = sin(x.Re());
		real cos_x = cos(x.Re());

		return dual2(sin_x,
			cos_x * x.Dual1(),
			-sin_x * square(x.Dual1()) + cos_x * x.Dual2());
	}


	/// Compute the cosine of a second order dual number
	inline dual2 cos(dual2 x) {

		real sin_x = sin(x.Re());
		real cos_x = cos(x.Re());

		return dual2(cos_x,
			-sin_x * x.Dual1(),
			-cos_x * square(x.Dual1()) - sin_x * x.Dual2());
	}


	/// Compute the tangent of a second order dual number
	inline dual2 tan(dual2 x) {

		real sin_x = sin(x.Re());
		real cos_x = cos(x.Re());

		if(cos_x == 0) {
			TH_MATH_ERROR("tan(dual2)", cos_x, DIV_BY_ZERO);
			return dual2(nan(), nan());
		}

		return dual2(tan(x.Re()),
			x.Dual1() / square(cos_x),
			2.0 * sin_x / cube(cos_x) * square(x.Dual1())
				+ x.Dual2() / square(cos_x));
	}


	/// Compute the cotangent of a second order dual number
	inline dual2 cot(dual2 x) {

		real sin_x = sin(x.Re());
		real cos_x = cos(x.Re());

		if(sin_x == 0) {
			TH_MATH_ERROR("cot(dual2)", sin_x, DIV_BY_ZERO);
			return dual2(nan(), nan());
		}

		return dual2(cot(x.Re()),
			-x.Dual1() / square(sin_x),
			-2.0 * cos_x / cube(sin_x) * square(x.Dual1())
				- x.Dual2() / square(sin_x));
	}


	/// Compute the exponential of a second order dual number
	inline dual2 exp(dual2 x) {
		real exp_x = exp(x.Re());
		return dual2(exp_x,
			x.Dual1() * exp_x,
			square(x.Dual1()) * exp_x + x.Dual2() * exp_x);
	}


	/// Compute the natural logarithm of a second order dual number
	inline dual2 ln(dual2 x) {

		if(x.Re() <= 0) {
			TH_MATH_ERROR("ln(dual2)", x.Re(), OUT_OF_DOMAIN);
			return dual2(nan(), nan());
		}

		return dual2(ln(x.Re()),
			x.Dual1() / x.Re(),
			-square(x.Dual1()) / square(x.Re()) + x.Dual2() / x.Re());
	}


	/// Compute the natural logarithm of a second order dual number
	inline dual2 log2(dual2 x) {

		if(x.Re() <= 0) {
			TH_MATH_ERROR("log2(dual2)", x.Re(), OUT_OF_DOMAIN);
			return dual2(nan(), nan());
		}

		return dual2(log2(x.Re()),
			x.Dual1() * LOG2E / x.Re(),
			-square(x.Dual1()) * LOG2E / square(x.Re())
				+ x.Dual2() * LOG2E / x.Re());
	}


	/// Compute the natural logarithm of a second order dual number
	inline dual2 log10(dual2 x) {

		if(x.Re() <= 0) {
			TH_MATH_ERROR("log10(dual2)", x.Re(), OUT_OF_DOMAIN);
			return dual2(nan(), nan());
		}

		return dual2(log10(x.Re()),
			x.Dual1() * LOG10E / x.Re(),
			-square(x.Dual1()) * LOG10E / square(x.Re())
				+ x.Dual2() * LOG10E / x.Re());
	}


	/// Compute the absolute value of a second order dual number
	inline dual2 abs(dual2 x) {
		return dual2(abs(x.Re()), x.Dual1() * sgn(x.Re()), 0);
	}


	/// Compute the arcsine of a second order dual number
	inline dual2 asin(dual2 x) {

		if(x.Re() >= 1) {
			TH_MATH_ERROR("asin(dual2)", x.Re(), OUT_OF_DOMAIN);
			return dual2(nan(), nan());
		}

		real sqrt_1mx2 = sqrt(1 - square(x.Re()));

		return dual2(asin(x.Re()),
			x.Dual1() / sqrt_1mx2,
			square(x.Dual1()) * x.Re() / ((1 - square(x.Re())) * sqrt_1mx2)
				+ x.Dual2() / sqrt_1mx2);
	} 


	/// Compute the arcosine of a second order dual number
	inline dual2 acos(dual2 x) {

		if(x.Re() >= 1) {
			TH_MATH_ERROR("acos(dual2)", x.Re(), OUT_OF_DOMAIN);
			return dual2(nan(), nan());
		}

		real sqrt_1mx2 = sqrt(1 - square(x.Re()));

		return dual2(acos(x.Re()),
			-x.Dual1() / sqrt_1mx2,
			-square(x.Dual1()) * x.Re() / ((1 - square(x.Re())) * sqrt_1mx2)
				- x.Dual2() / sqrt_1mx2);
	}


	/// Compute the arctangent of a second order dual number
	inline dual2 atan(dual2 x) {
		return dual2(atan(x.Re()),
			x.Dual1() / (1 + square(x.Re())),
			square(x.Dual1()) * 2 * x.Re() / square(1 + square(x.Re()))
				+ x.Dual2() / (1 + square(x.Re())));
	}

}


#endif
