
///
/// @file Functions defined on multidual numbers for
/// automatic differentiation of multivariable real functions


#ifndef UROBORO_MULTIDUAL_FUNCTIONS_H
#define UROBORO_MULTIDUAL_FUNCTIONS_H

#include "./multidual.h"
#include "../real_analysis.h"


namespace uroboro {


	/// Return the square of a multidual number
	template<unsigned int N>
	multidual<N> square(multidual<N> x) {
		return x * x;
	}


	/// Return the cube of a multidual number
	template<unsigned int N>
	multidual<N> cube(multidual<N> x) {
		return x * x * x;
	}


	/// Compute the n-th power of a multidual number
	template<unsigned int N>
	multidual<N> pow(multidual<N> x, int n) {
		real pow_n_1_x = pow(x.Re(), n - 1);
		return multidual<N>(pow_n_1_x * x.Re(), x.Dual() * pow_n_1_x * n);
	}


	/// Compute the square root of a multidual number
	template<unsigned int N>
	multidual<N> sqrt(multidual<N> x) {
		real sqrt_x = sqrt(x.Re());
		return multidual<N>(sqrt_x, x.Dual() * 0.5 / sqrt_x);
	}


	/// Compute the sine of a multidual number
	template<unsigned int N>
	multidual<N> sin(multidual<N> x) {
		return multidual<N>(sin(x.Re()), x.Dual() * cos(x.Re()));
	}


	/// Compute the cosine of a multidual number
	template<unsigned int N>
	multidual<N> cos(multidual<N> x) {
		return multidual<N>(cos(x.Re()), x.Dual() * -sin(x.Re()));
	}


	/// Compute the tangent of a multidual number
	template<unsigned int N>
	multidual<N> tan(multidual<N> x) {
		return multidual<N>(tan(x.Re()), x.Dual() / square(cos(x.Re())));
	}


	/// Compute the cotangent of a multidual number
	template<unsigned int N>
	multidual<N> cot(multidual<N> x) {
		return multidual<N>(cot(x.Re()), x.Dual() * (-1 / square(sin(x.Re()))));
	}


	/// Compute the exponential of a multidual number
	template<unsigned int N>
	multidual<N> exp(multidual<N> x) {
		return multidual<N>(exp(x.Re()), x.Dual() * exp(x.Re()));
	}


	/// Compute the natural logarithm of a multidual number
	template<unsigned int N>
	multidual<N> ln(multidual<N> x) {
		return multidual<N>(ln(x.Re()), x.Dual() / x.Re());
	}


	/// Compute the natural logarithm of a multidual number
	template<unsigned int N>
	multidual<N> log2(multidual<N> x) {
		return multidual<N>(log2(x.Re()), x.Dual() * LOG2E / x.Re());
	}


	/// Compute the natural logarithm of a multidual number
	template<unsigned int N>
	multidual<N> log10(multidual<N> x) {
		return multidual<N>(log10(x.Re()), x.Dual() * LOG10E / x.Re());
	}


	/// Compute the absolute value of a multidual number
	template<unsigned int N>
	multidual<N> abs(multidual<N> x) {
		return multidual<N>(abs(x.Re()), x.Dual() * sgn(x.Re()));
	}


	/// Compute the arcsine of a multidual number
	template<unsigned int N>
	multidual<N> asin(multidual<N> x) {
		return multidual<N>(asin(x.Re()), x.Dual() / sqrt(1 - square(x.Re())));
	}


	/// Compute the arcsine of a multidual number
	template<unsigned int N>
	multidual<N> acos(multidual<N> x) {
		return multidual<N>(acos(x.Re()), x.Dual() * (-1 / sqrt(1 - square(x.Re()))));
	}

	/// Compute the arcsine of a multidual number
	template<unsigned int N>
	multidual<N> atan(multidual<N> x) {
		return multidual<N>(atan(x.Re()), x.Dual() / (1 + square(x.Re())));
	}

}


#endif
