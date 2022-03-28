
/// 
/// @file extrema.h Extrema approximation of real functions
/// 

#ifndef UROBORO_EXTREMA_H
#define UROBORO_EXTREMA_H

#include "../constants.h"
#include "./roots.h"


namespace uroboro {


	/// Approximate a function maximum given the function and the first two derivatives using Newton-Raphson
	inline real approx_max_newton(
		real_function f, real_function Df, real_function D2f,
		real guess = 0) {

		real z = approx_root_newton(Df, D2f, guess);

		if(D2f(z) > 0) {
			UMATH_ERROR("approx_max_newton", z, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return z;
	}


	/// Approximate a function minimum given the function and the first two derivatives using Newton-Raphson
	inline real approx_min_newton(
		real_function f, real_function Df,
		real_function D2f, real guess = 0) {

		real z = approx_root_newton(Df, D2f, guess);

		if(D2f(z) < 0) {
			UMATH_ERROR("approx_min_newton", z, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return z;
	}


	/// Approximate a function maximum inside an interval given
	/// the function and its first derivative using bisection on the derivative
	inline real approx_max_bisection(
		real_function f, real_function Df,
		real a, real b) {

		real z = approx_root_bisection(Df, a, b);

		if(approx_derivative(Df, z) > 0) {
			UMATH_ERROR("approx_max_bisection", z, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return z;
	}


	/// Approximate a function minimum inside an interval given the function
	/// and its first derivative using bisection on the derivative
	inline real approx_min_bisection(real_function f, real_function Df, real a, real b) {

		real z = approx_root_bisection(Df, a, b);

		if(approx_derivative(Df, z) < 0) {
			UMATH_ERROR("approx_min_bisection", z, NO_ALGO_CONVERGENCE);
			return z;
		}

		return z;
	}

}

#endif
