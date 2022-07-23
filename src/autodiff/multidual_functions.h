
///
/// @file multidual_functions.h Functions defined on multidual numbers for
/// automatic differentiation of multivariable real functions


#ifndef THEORETICA_MULTIDUAL_FUNCTIONS_H
#define THEORETICA_MULTIDUAL_FUNCTIONS_H

#include "./multidual.h"
#include "../core/real_analysis.h"


namespace theoretica {


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

		if(sqrt_x == 0) {
			TH_MATH_ERROR("sqrt(multidual)", sqrt_x, DIV_BY_ZERO);
			return multidual<N>(nan(), vec<N>(nan()));
		}

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

		real cos_x = cos(x.Re());

		if(cos_x == 0) {
			TH_MATH_ERROR("tan(multidual)", cos_x, DIV_BY_ZERO);
			return multidual<N>(nan(), vec<N>(nan()));
		}

		return multidual<N>(tan(x.Re()), x.Dual() / square(cos_x));
	}


	/// Compute the cotangent of a multidual number
	template<unsigned int N>
	multidual<N> cot(multidual<N> x) {

		real sin_x = sin(x.Re());

		if(sin_x == 0) {
			TH_MATH_ERROR("cot(multidual)", sin_x, DIV_BY_ZERO);
			return multidual<N>(nan(), vec<N>(nan()));
		}

		return multidual<N>(cot(x.Re()), x.Dual() * (-1 / square(sin_x)));
	}


	/// Compute the exponential of a multidual number
	template<unsigned int N>
	multidual<N> exp(multidual<N> x) {
		real exp_x = exp(x.Re());
		return multidual<N>(exp_x, x.Dual() * exp_x);
	}


	/// Compute the natural logarithm of a multidual number
	template<unsigned int N>
	multidual<N> ln(multidual<N> x) {

		if(x.Re() <= 0) {
			TH_MATH_ERROR("ln(multidual)", x.Re(), OUT_OF_DOMAIN);
			return multidual<N>(nan(), vec<N>(nan()));
		}

		return multidual<N>(ln(x.Re()), x.Dual() / x.Re());
	}


	/// Compute the natural logarithm of a multidual number
	template<unsigned int N>
	multidual<N> log2(multidual<N> x) {

		if(x.Re() <= 0) {
			TH_MATH_ERROR("log2(multidual)", x.Re(), OUT_OF_DOMAIN);
			return multidual<N>(nan(), vec<N>(nan()));
		}

		return multidual<N>(log2(x.Re()), x.Dual() * LOG2E / x.Re());
	}


	/// Compute the natural logarithm of a multidual number
	template<unsigned int N>
	multidual<N> log10(multidual<N> x) {

		if(x.Re() <= 0) {
			TH_MATH_ERROR("log10(multidual)", x.Re(), OUT_OF_DOMAIN);
			return multidual<N>(nan(), vec<N>(nan()));
		}

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

		if(x.Re() >= 1) {
			TH_MATH_ERROR("asin(multidual)", x.Re(), OUT_OF_DOMAIN);
			return multidual<N>(nan(), vec<N>(nan()));
		}

		return multidual<N>(asin(x.Re()), x.Dual() / sqrt(1 - square(x.Re())));
	}


	/// Compute the arccosine of a multidual number
	template<unsigned int N>
	multidual<N> acos(multidual<N> x) {

		if(x.Re() >= 1) {
			TH_MATH_ERROR("acos(multidual)", x.Re(), OUT_OF_DOMAIN);
			return multidual<N>(nan(), vec<N>(nan()));
		}

		return multidual<N>(acos(x.Re()), x.Dual() * (-1 / sqrt(1 - square(x.Re()))));
	}

	/// Compute the arctangent of a multidual number
	template<unsigned int N>
	multidual<N> atan(multidual<N> x) {
		return multidual<N>(atan(x.Re()), x.Dual() / (1 + square(x.Re())));
	}


	/// Compute the hyperbolic sine of a multidual number
	template<unsigned int N>
	multidual<N> sinh(multidual<N> x) {

		real exp_x = exp(x.Re());
		return multidual<N>((exp_x - 1.0 / exp_x) / 2.0, x.Dual() * (exp_x + 1.0 / exp_x) / 2.0);
	}


	/// Compute the hyperbolic cosine of a multidual number
	template<unsigned int N>
	multidual<N> cosh(multidual<N> x) {

		real exp_x = exp(x.Re());
		return multidual<N>((exp_x + 1.0 / exp_x) / 2.0, x.Dual() * (exp_x - 1.0 / exp_x) / 2.0);
	}


	/// Compute the hyperbolic tangent of a multidual number
	template<unsigned int N>
	multidual<N> tanh(multidual<N> x) {

		real exp_x = exp(x.Re());
		return multidual<N>(
			(exp_x - 1.0 / exp_x) / (exp_x + 1.0 / exp_x),
			x.Dual() / square(exp_x + 1.0 / exp_x));
	}

}


#endif
