
/// 
/// @file extrema.h Extrema approximation of real functions
/// 

#ifndef THEORETICA_EXTREMA_H
#define THEORETICA_EXTREMA_H

#include "../core/constants.h"
#include "./roots.h"


namespace theoretica {


	/// Approximate a function maximum using the Golden Section search algorithm
	inline real approx_max_goldensection(real_function f, real a, real b) {

		if(a > b) {
			TH_MATH_ERROR("approx_max_goldensection", b, INVALID_ARGUMENT);
			return nan();
		}

		real x1 = a;
		real x2 = b;
		real x3 = b - (b - a) / PHI;
		real x4 = a + (b - a) / PHI;

		unsigned int iter = 0;

		while(abs(x2 - x1) > ROOT_APPROX_TOL && iter <= MAX_GOLDENSECTION_ITER) {

			if(f(x3) > f(x4)) {
				x2 = x4;
			} else {
				x1 = x3;
			}

			x3 = x2 - (x2 - x1) / PHI;
			x4 = x1 + (x2 - x1) / PHI;

			iter++;
		}

		if(iter > MAX_GOLDENSECTION_ITER) {
			TH_MATH_ERROR("approx_max_goldensection", iter, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return (x2 + x1) / 2.0;
	}


	/// Approximate a function minimum using the Golden Section search algorithm
	inline real approx_min_goldensection(real_function f, real a, real b) {

		if(a > b) {
			TH_MATH_ERROR("approx_min_goldensection", b, INVALID_ARGUMENT);
			return nan();
		}

		real x1 = a;
		real x2 = b;
		real x3 = b - (b - a) / PHI;
		real x4 = a + (b - a) / PHI;

		unsigned int iter = 0;

		while(abs(x2 - x1) > ROOT_APPROX_TOL && iter <= MAX_GOLDENSECTION_ITER) {

			if(f(x3) < f(x4)) {
				x2 = x4;
			} else {
				x1 = x3;
			}

			x3 = x2 - (x2 - x1) / PHI;
			x4 = x1 + (x2 - x1) / PHI;

			iter++;
		}

		if(iter > MAX_GOLDENSECTION_ITER) {
			TH_MATH_ERROR("approx_min_goldensection", iter, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return (x2 + x1) / 2.0;
	}


	/// Approximate a function maximum given the function and the first two derivatives using Newton-Raphson
	inline real approx_max_newton(
		real_function f, real_function Df, real_function D2f,
		real guess = 0) {

		real z = approx_root_newton(Df, D2f, guess);

		if(D2f(z) > 0) {
			TH_MATH_ERROR("approx_max_newton", z, NO_ALGO_CONVERGENCE);
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
			TH_MATH_ERROR("approx_min_newton", z, NO_ALGO_CONVERGENCE);
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

		if(deriv_central(Df, z) > 0) {
			TH_MATH_ERROR("approx_max_bisection", z, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return z;
	}


	/// Approximate a function minimum inside an interval given the function
	/// and its first derivative using bisection on the derivative
	inline real approx_min_bisection(real_function f, real_function Df, real a, real b) {

		real z = approx_root_bisection(Df, a, b);

		if(deriv_central(Df, z) < 0) {
			TH_MATH_ERROR("approx_min_bisection", z, NO_ALGO_CONVERGENCE);
			return z;
		}

		return z;
	}

}

#endif
