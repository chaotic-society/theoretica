
///
/// @file dual_functions.h Functions defined on dual numbers for
/// automatic differentiation.
///
/// Dual numbers can be easily used for
/// automatic differentiation, as they behave
/// under addition, multiplication and division
/// as the first derivative.
///
/// Normal operations can be performed and the result
/// will have a real part equal to the function evaluated
/// for the given argument and a "dual" part equal
/// to the first derivative evaluated for the given argument.


#ifndef THEORETICA_DUAL_FUNCTIONS_H
#define THEORETICA_DUAL_FUNCTIONS_H

#include "./dual.h"
#include "../core/real_analysis.h"


namespace theoretica {


	/// Return the square of a dual number
	inline dual square(dual x) {
		return x * x;
	}


	/// Return the cube of a dual number
	inline dual cube(dual x) {
		return x * x * x;
	}


	/// Compute the n-th power of a dual number
	inline dual pow(dual x, int n) {
		real pow_n_1_x = pow(x.Re(), n - 1);
		return dual(pow_n_1_x * x.Re(), pow_n_1_x * n * x.Dual());
	}


	/// Compute the square root of a dual number
	inline dual sqrt(dual x) {

		real sqrt_x = sqrt(x.Re());

		if(sqrt_x == 0) {
			TH_MATH_ERROR("sqrt(dual)", sqrt_x, DIV_BY_ZERO);
			return dual(nan(), nan());
		}

		return dual(sqrt_x, 0.5 / sqrt_x * x.Dual());
	}


	/// Compute the sine of a dual number
	inline dual sin(dual x) {
		return dual(sin(x.Re()), cos(x.Re()) * x.Dual());
	}


	/// Compute the cosine of a dual number
	inline dual cos(dual x) {
		return dual(cos(x.Re()), -sin(x.Re()) * x.Dual());
	}


	/// Compute the tangent of a dual number
	inline dual tan(dual x) {

		real cos_x = cos(x.Re());

		if(cos_x == 0) {
			TH_MATH_ERROR("tan(dual)", cos_x, DIV_BY_ZERO);
			return dual(nan(), nan());
		}

		return dual(tan(x.Re()), x.Dual() / square(cos_x));
	}


	/// Compute the cotangent of a dual number
	inline dual cot(dual x) {

		real sin_x = sin(x.Re());

		if(sin_x == 0) {
			TH_MATH_ERROR("cot(dual)", sin_x, DIV_BY_ZERO);
			return dual(nan(), nan());
		}

		return dual(cot(x.Re()), -x.Dual() / square(sin_x));
	}


	/// Compute the exponential of a dual number
	inline dual exp(dual x) {
		real exp_x = exp(x.Re());
		return dual(exp_x, x.Dual() * exp_x);
	}


	/// Compute the natural logarithm of a dual number
	inline dual ln(dual x) {

		if(x.Re() <= 0) {
			TH_MATH_ERROR("ln(dual)", x.Re(), OUT_OF_DOMAIN);
			return dual(nan(), nan());
		}

		return dual(ln(x.Re()), x.Dual() / x.Re());
	}


	/// Compute the natural logarithm of a dual number
	inline dual log2(dual x) {

		if(x.Re() <= 0) {
			TH_MATH_ERROR("log2(dual)", x.Re(), OUT_OF_DOMAIN);
			return dual(nan(), nan());
		}

		return dual(log2(x.Re()), x.Dual() * LOG2E / x.Re());
	}


	/// Compute the natural logarithm of a dual number
	inline dual log10(dual x) {

		if(x.Re() <= 0) {
			TH_MATH_ERROR("log10(dual)", x.Re(), OUT_OF_DOMAIN);
			return dual(nan(), nan());
		}

		return dual(log10(x.Re()), x.Dual() * LOG10E / x.Re());
	}


	/// Compute the absolute value of a dual number
	inline dual abs(dual x) {
		return dual(abs(x.Re()), x.Dual() * sgn(x.Re()));
	}


	/// Compute the arcsine of a dual number
	inline dual asin(dual x) {

		if(x.Re() >= 1) {
			TH_MATH_ERROR("asin(dual)", x.Re(), OUT_OF_DOMAIN);
			return dual(nan(), nan());
		}

		return dual(asin(x.Re()), x.Dual() / sqrt(1 - square(x.Re())));
	}


	/// Compute the arccosine of a dual number
	inline dual acos(dual x) {

		if(x.Re() >= 1) {
			TH_MATH_ERROR("acos(dual)", x.Re(), OUT_OF_DOMAIN);
			return dual(nan(), nan());
		}

		return dual(acos(x.Re()), -x.Dual() / sqrt(1 - square(x.Re())));
	}


	/// Compute the arctangent of a dual number
	inline dual atan(dual x) {
		return dual(atan(x.Re()), x.Dual() / (1 + square(x.Re())));
	}


	/// Compute the hyperbolic sine of a dual number
	inline dual sinh(dual x) {

		real exp_x = exp(x.Re());
		return dual((exp_x - 1.0 / exp_x) / 2.0, x.Dual() * (exp_x + 1.0 / exp_x) / 2.0);
	}


	/// Compute the hyperbolic cosine of a dual number
	inline dual cosh(dual x) {

		real exp_x = exp(x.Re());
		return dual((exp_x + 1.0 / exp_x) / 2.0, x.Dual() * (exp_x - 1.0 / exp_x) / 2.0);
	}


	/// Compute the hyperbolic tangent of a dual number
	inline dual tanh(dual x) {

		real exp_x = exp(x.Re());
		return dual(
			(exp_x - 1.0 / exp_x) / (exp_x + 1.0 / exp_x),
			x.Dual() / square(exp_x + 1.0 / exp_x));
	}

}


#endif
