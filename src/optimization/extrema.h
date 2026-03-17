
/// 
/// @file extrema.h Extrema approximation of real functions
/// 

#ifndef THEORETICA_EXTREMA_H
#define THEORETICA_EXTREMA_H

#include "../core/constants.h"
#include "./roots.h"


namespace theoretica {


	/// Approximate a function maximum using the
	/// Golden Section search algorithm.
	///
	/// @param f The function to search a local maximum of.
	/// @param a The lower extreme of the search interval.
	/// @param b The upper extreme of the search interval.
	/// @return The coordinate of the local maximum.
	template<typename RealFunction>
	inline iter_result<real> maximize_goldensection(
		RealFunction f, real a, real b,
		real tolerance = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_GOLDENSECTION_ITER) {

		if(a > b) {
			TH_MATH_ERROR("maximize_goldensection", b, MathError::InvalidArgument);
			return iter_result<real>(ConvergenceStatus::InvalidInput);
		}

		real x1 = a;
		real x2 = b;
		real x3 = b - (b - a) / PHI;
		real x4 = a + (b - a) / PHI;

		unsigned int iter = 0;

		// Keep iterating until the bracketing interval is closer than the tolerance,
		// or until the maximum number of iterations is reached.
		while(abs(x2 - x1) > tolerance && iter <= max_iter) {

			if(f(x3) > f(x4)) {
				x2 = x4;
			} else {
				x1 = x3;
			}

			x3 = x2 - (x2 - x1) / PHI;
			x4 = x1 + (x2 - x1) / PHI;

			iter++;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("maximize_goldensection", iter, MathError::NoConvergence);
			return iter_result<real>(ConvergenceStatus::MaxIterations, iter, abs(x2 - x1) / 2.0);
		}

		return iter_result<real>((x2 + x1) / 2.0, iter, abs(x2 - x1) / 2.0);
	}


	/// Approximate a function minimum using the
	/// Golden Section search algorithm.
	///
	/// @param f The function to search a local minimum of.
	/// @param a The lower extreme of the search interval.
	/// @param b The upper extreme of the search interval.
	/// @return The coordinate of the local minimum.
	template<typename RealFunction>
	inline iter_result<real> minimize_goldensection(
		RealFunction f, real a, real b,
		real tolerance = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_GOLDENSECTION_ITER) {

		if(a > b) {
			TH_MATH_ERROR("minimize_goldensection", b, MathError::InvalidArgument);
			return iter_result<real>(ConvergenceStatus::InvalidInput);
		}

		real x1 = a;
		real x2 = b;
		real x3 = b - (b - a) / PHI;
		real x4 = a + (b - a) / PHI;

		unsigned int iter = 0;

		while(abs(x2 - x1) > tolerance && iter <= max_iter) {

			if(f(x3) < f(x4)) {
				x2 = x4;
			} else {
				x1 = x3;
			}

			x3 = x2 - (x2 - x1) / PHI;
			x4 = x1 + (x2 - x1) / PHI;

			iter++;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("minimize_goldensection", iter, MathError::NoConvergence);
			return iter_result<real>(ConvergenceStatus::MaxIterations, iter, abs(x2 - x1) / 2.0);
		}

		return iter_result<real>((x2 + x1) / 2.0, iter, abs(x2 - x1) / 2.0);
	}


	/// Approximate a function maximum given the function and the
	/// first two derivatives using Newton-Raphson's method
	/// to find a root of the derivative.
	///
	/// @param f The function to search a local maximum of.
	/// @param Df The derivative of the function.
	/// @param D2f The second derivative of the function.
	/// @param guess The initial guess (defaults to 0).
	/// @return The coordinate of the local maximum.
	template<typename RealFunction>
	inline iter_result<real> maximize_newton(
		RealFunction f, RealFunction Df, RealFunction D2f,
		real guess = 0.0, real tolerance = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_NEWTON_ITER) {

		iter_result<real> z = root_newton(Df, D2f, guess, tolerance, max_iter);

		if(D2f(z) > 0) {
			TH_MATH_ERROR("maximize_newton", z, MathError::NoConvergence);
			return iter_result<real>(z.status, z.iterations, z.residual);
		}

		return iter_result<real>(z, z.iterations, z.residual);
	}


	/// Approximate a function minimum given the function
	/// and the first two derivatives using Newton-Raphson's
	/// method to find a root of the derivative.
	///
	/// @param f The function to search a local minimum of.
	/// @param Df The derivative of the function.
	/// @param D2f The second derivative of the function.
	/// @param guess The initial guess (defaults to 0).
	/// @return The coordinate of the local minimum.
	template<typename RealFunction>
	inline iter_result<real> minimize_newton(
		RealFunction f, RealFunction Df,
		RealFunction D2f, real guess = 0,
		real tolerance = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_NEWTON_ITER) {

		iter_result<real> z = root_newton(Df, D2f, guess, tolerance, max_iter);

		if(D2f(z) < 0) {
			TH_MATH_ERROR("minimize_newton", z, MathError::NoConvergence);
			return iter_result<real>(z.status, z.iterations, z.residual);
		}

		return iter_result<real>(z, z.iterations, z.residual);
	}


	/// Approximate a function maximum inside an interval given
	/// the function and its first derivative using bisection
	/// on the derivative.
	///
	/// @param f The function to search a local minimum of.
	/// @param Df The derivative of the function.
	/// @param a The lower extreme of the search interval.
	/// @param b The upper extreme of the search interval.
	/// @return The coordinate of the local maximum.
	template<typename RealFunction>
	inline iter_result<real> maximize_bisection(
		RealFunction f, RealFunction Df,
		real a, real b,
		real tolerance = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_BISECTION_ITER) {

		iter_result<real> z = root_bisect(Df, a, b, tolerance, max_iter);

		if(deriv_central(Df, z) > 0) {
			TH_MATH_ERROR("maximize_bisection", z, MathError::NoConvergence);
			return iter_result<real>(z.status, z.iterations, z.residual);
		}

		return iter_result<real>(z, z.iterations, z.residual);
	}


	/// Approximate a function minimum inside an interval
	/// given the function and its first derivative using
	/// bisection on the derivative.
	///
	/// @param f The function to search a local minimum of.
	/// @param Df The derivative of the function.
	/// @param a The lower extreme of the search interval.
	/// @param b The upper extreme of the search interval.
	/// @return The coordinate of the local minimum.
	template<typename RealFunction>
	inline iter_result<real> minimize_bisection(
		RealFunction f, RealFunction Df, real a, real b,
		real tolerance = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_BISECTION_ITER) {

		iter_result<real> z = root_bisect(Df, a, b, tolerance, max_iter);

		if(deriv_central(Df, z) < 0) {
			TH_MATH_ERROR("minimize_bisection", z, MathError::NoConvergence);
			return iter_result<real>(z.status, z.iterations, z.residual);
		}

		return iter_result<real>(z, z.iterations, z.residual);
	}

}

#endif
