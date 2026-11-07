
///
/// @file multi_roots.h Numerical methods for multivariate root finding
///

#ifndef THEORETICA_MULTI_ROOTS_H
#define THEORETICA_MULTI_ROOTS_H

#include "../autodiff/autodiff.h"
#include "../core/iter_result.h"


namespace theoretica {


	/// Approximate the root of a multivariate function
	/// using Newton's method with pure Jacobian.
	///
	/// @note For automatic differentiation, only vec<> types are supported
	/// for argument and return vectors in f (either fixed or dynamic size).
	///
	/// @param f The function to find the root of
	/// @param guess The first guess (defaults to the origin)
	/// @param tolerance The tolerance over the final result
	/// (defaults to OPTIMIZATION_MINGRAD_TOLERANCE)
	/// @param max_iter The maximum number of iterations before
	/// stopping the algorithm
	/// @result The computed vector at which f is approximately zero
#ifdef THEORETICA_HAS_CPP20
	template <
		VectorType Vector,
		VectorType ReturnVector = Vector,
		autodiff::ADVectorField VectorField
	>
#else
	template <
		typename Vector = vec<real>,
		typename ReturnVector = Vector,
		typename VectorField,
		autodiff::enable_vector_field<VectorField> = true
	>
#endif
	inline iter_result<ReturnVector> multiroot_newton(
		VectorField f,
		Vector guess,
		real tolerance = OPTIMIZATION_MINGRAD_TOLERANCE,
		unsigned int max_iter = OPTIMIZATION_MINGRAD_ITER
	) {

		// Current position
		ReturnVector x = guess;

		// Number of iterations
		unsigned int iter = 0;

		// Extract the size of the vector type
		constexpr size_t N = ReturnVector::size_argument;

		// The current value of f(x)
		ReturnVector f_x = multidual<N>::extract_real(
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
			x -= algebra::solve(J, f_x);
			iter++;
		}

		if(iter > max_iter) {

			TH_MATH_ERROR("multiroot_newton", iter, MathError::NoConvergence);
			algebra::vec_error(x);
			return iter_result<ReturnVector>(
				x, ConvergenceStatus::MaxIterations, iter, algebra::norm(f_x)
			);
		}

		return iter_result<ReturnVector>(x, iter, algebra::norm(f_x));
	}

}

#endif
