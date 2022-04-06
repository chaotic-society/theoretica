
///
/// @file Functions defined on dual numbers for
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


#ifndef UROBORO_DUAL_FUNCTIONS_H
#define UROBORO_DUAL_FUNCTIONS_H

#include "./dual.h"
#include "../real_analysis.h"


namespace uroboro {


	/// Return the square of a dual number
	dual square(dual x) {
		return x * x;
	}


	/// Return the cube of a dual number
	dual cube(dual x) {
		return x * x * x;
	}


	/// Compute the n-th power of a dual number
	dual pow(dual x, int n) {
		real pow_n_1_x = pow(x.Re(), n - 1);
		return dual(pow_n_1_x * x.Re(), pow_n_1_x * n * x.Dual());
	}


	/// Compute the square root of a dual number
	dual sqrt(dual x) {
		real sqrt_x = sqrt(x.Re());
		return dual(sqrt_x, 0.5 / sqrt_x * x.Dual());
	}


	/// Compute the sine of a dual number
	dual sin(dual x) {
		return dual(sin(x.Re()), cos(x.Re()) * x.Dual());
	}


	/// Compute the cosine of a dual number
	dual cos(dual x) {
		return dual(cos(x.Re()), -sin(x.Re()) * x.Dual());
	}


	/// Compute the tangent of a dual number
	dual tan(dual x) {
		return dual(tan(x.Re()), x.Dual() / square(cos(x.Re())));
	}


	/// Compute the cotangent of a dual number
	dual cot(dual x) {
		return dual(cot(x.Re()), -x.Dual() / square(sin(x.Re())));
	}


	/// Compute the exponential of a dual number
	dual exp(dual x) {
		return dual(exp(x.Re()), x.Dual() * exp(x.Re()));
	}


	/// Compute the natural logarithm of a dual number
	dual ln(dual x) {
		return dual(ln(x.Re()), x.Dual() / x.Re());
	}


	/// Compute the natural logarithm of a dual number
	dual log2(dual x) {
		return dual(log2(x.Re()), x.Dual() * LOG2E / x.Re());
	}


	/// Compute the natural logarithm of a dual number
	dual log10(dual x) {
		return dual(log10(x.Re()), x.Dual() * LOG10E / x.Re());
	}


	/// Compute the absolute value of a dual number
	dual abs(dual x) {
		return dual(abs(x.Re()), x.Dual() * sgn(x.Re()));
	}


	/// Compute the arcsine of a dual number
	dual asin(dual x) {
		return dual(asin(x.Re()), x.Dual() / sqrt(1 - square(x.Re())));
	}


	/// Compute the arcsine of a dual number
	dual acos(dual x) {
		return dual(acos(x.Re()), -x.Dual() / sqrt(1 - square(x.Re())));
	}


	/// Compute the arcsine of a dual number
	dual atan(dual x) {
		return dual(atan(x.Re()), x.Dual() / (1 + square(x.Re())));
	}

}


#endif
