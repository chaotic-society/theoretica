
///
/// @file multi_roots.h Numerical methods for multivariate root finding
///

#ifndef THEORETICA_MULTI_ROOTS_H
#define THEORETICA_MULTI_ROOTS_H

#include "../autodiff/autodiff.h"


namespace theoretica {


	/// Approximate the root of a multivariate function
	/// using Newton's method with pure Jacobian.
	///
	/// @param f The function to find the root of
	/// @param guess The first guess (defaults to the origin)
	/// @param tolerance The tolerance over the final result
	/// (defaults to MINGRAD_TOLERANCE)
	/// @param max_iter The maximum number of iterations before
	/// stopping the algorithm
	/// @result The computed vector at which f is approximately zero
	template<unsigned int N>
	inline vec<real, N> multiroot_newton(
		d_vec<N>(*f)(d_vec<N>),
		vec<real, N> guess = vec<real, N>(0),
		real tolerance = MINGRAD_TOLERANCE,
		unsigned int max_iter = MINGRAD_MAX_ITER) {

		// Current position
		vec<real, N> x = guess;

		// Number of iterations
		unsigned int iter = 0;

		// The current value of f(x)
		vec<real, N> f_x = multidual<N>::extract_real(
			f(multidual<N>::make_argument(x))
		);

		// Jacobian matrix
		mat<real, N, N> J;

		while(f_x.sqr_norm() > square(tolerance) && iter <= max_iter) {

			// Compute the function value and Jacobian
			// using automatic differentiation
			multidual<N>::extract(
				f(multidual<N>::make_argument(x)), f_x, J
			);

			// Update the current best guess
			x = x - J.inverse() * f_x;
			iter++;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("multi_root_newton", iter, NO_ALGO_CONVERGENCE);
			return vec<real, N>(nan());
		}

		return x;
	}

}

#endif
