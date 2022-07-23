
///
/// @file multi_extrema.h Search of extrema of multivariate functions
///

#ifndef THEORETICA_MULTI_EXTREMA_H
#define THEORETICA_MULTI_EXTREMA_H

#include "../core/constants.h"
#include "../autodiff/autodiff.h"
#include "./extrema.h"


namespace theoretica {


	/// Find a local minimum of the given multivariate function
	/// using fixed-step gradient descent
	///
	/// @param f The function to minimize
	/// @param guess The initial guess, defaults to {0, 0}
	/// @param gamma The fixed step size, defaults to MINGRAD_GAMMA
	/// @param tolerance The maximum magnitude of the gradient to stop
	/// the algorithm at, defaults to MINGRAD_TOLERANCE.
	/// @param max_iter The maximum number of iterations to perform before
	/// stopping execution of the routine.
	/// @return The coordinates of the local minimum
	template<unsigned int N>
	inline vec<N> minimize_grad(
		multidual<N>(*f)(vec<N, multidual<N>>),
		vec<N> guess = vec<N>(0),
		real gamma = MINGRAD_GAMMA,
		real tolerance = MINGRAD_TOLERANCE,
		unsigned int max_iter = MINGRAD_MAX_ITER) {

		if(gamma >= 0) {
			TH_MATH_ERROR("minimize_grad", gamma, INVALID_ARGUMENT);
			return vec<N>(nan());
		}

		vec<N> x = guess;
		vec<N> grad;
		unsigned int iter = 0;

		do {

			grad = gradient(f, x);
			x += gamma * grad;
			iter++;

		} while(grad.lenght() > MINGRAD_TOLERANCE && iter <= max_iter);

		if(iter > max_iter) {
			TH_MATH_ERROR("minimize_grad", iter, NO_ALGO_CONVERGENCE);
			return vec<N>(nan());
		}

		return x;
	}


	/// Find a local maximum of the given multivariate function
	/// using fixed-step gradient descent
	///
	/// @param f The function to maximize
	/// @param guess The initial guess, defaults to {0, 0}
	/// @param gamma The fixed step size, defaults to MINGRAD_GAMMA
	/// @param tolerance The maximum magnitude of the gradient to stop
	/// the algorithm at, defaults to MINGRAD_TOLERANCE.
	/// @param max_iter The maximum number of iterations to perform before
	/// stopping execution of the routine.
	/// @return The coordinates of the local maximum
	template<unsigned int N>
	inline vec<N> maximize_grad(
		multidual<N>(*f)(vec<N, multidual<N>>),
		vec<N> guess = vec<N>(0),
		real gamma = MINGRAD_GAMMA,
		real tolerance = MINGRAD_TOLERANCE,
		unsigned int max_iter = MINGRAD_MAX_ITER) {

		return minimize_grad(f, guess, -gamma, tolerance, max_iter);
	}


	/// Find a local minimum of the given multivariate function
	/// using gradient descent with linear search
	///
	/// @param f The function to minimize
	/// @param guess The initial guess, defaults to {0, 0}
	/// @param tolerance The maximum magnitude of the gradient to stop
	/// the algorithm at, defaults to MINGRAD_TOLERANCE.
	/// @param max_iter The maximum number of iterations to perform before
	/// stopping execution of the routine.
	/// @return The coordinates of the local minimum
	template<unsigned int N>
	inline vec<N> minimize_lingrad(
		multidual<N>(*f)(vec<N, multidual<N>>),
		vec<N> guess = vec<N>(0),
		real tolerance = MINGRAD_TOLERANCE,
		unsigned int max_iter = MINGRAD_MAX_ITER) {

		vec<N> x = guess;
		vec<N> grad;
		unsigned int iter = 0;

		do {

			grad = gradient(f, x);

			// Minimize f(x + gamma * gradient) in [-1, 0]
			// using Golden Section extrema search
			real gamma = approx_min_goldensection(
				[f, x, grad](real gamma){
					return f(
						multidual<N>::pack_function_arg(x)
						+ multidual<N>::pack_function_arg(grad) * gamma).Re();
				},
				-1, 0);

			// Fallback value
			if(-gamma <= MACH_EPSILON)
				gamma = MINGRAD_GAMMA;

			x += gamma * grad;
			iter++;

		} while(grad.lenght() > MINGRAD_TOLERANCE && iter <= max_iter);

		if(iter > max_iter) {
			TH_MATH_ERROR("minimize_lingrad", iter, NO_ALGO_CONVERGENCE);
			return vec<N>(nan());
		}

		return x;
	}


	/// Use the best available algorithm to find a local
	/// minimum of the given multivariate function
	///
	/// @param f The function to minimize
	/// @param guess The initial guess
	/// @return The coordinates of the local minimum, 
	/// NaN if the algorithm did not converge.
	template<unsigned int N>
	inline vec<N> minimize(
		multidual<N>(*f)(vec<N, multidual<N>>),
		vec<N> guess = vec<N>(0)) {

		return minimize_lingrad(f, guess);
	}

}


#endif
