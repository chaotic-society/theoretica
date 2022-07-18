
///
/// @file multi_extrema.h Search of extrema of multivariate functions
///

#ifndef THEORETICA_MULTI_EXTREMA_H
#define THEORETICA_MULTI_EXTREMA_H

#include "../core/constants.h"
#include "../autodiff/autodiff.h"


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
	inline vec2 minimize_grad(
		multidual<2>(*f)(vec<2, multidual<2>>),
		vec2 guess = {0, 0},
		real gamma = MINGRAD_GAMMA,
		real tolerance = MINGRAD_TOLERANCE,
		unsigned int max_iter = MINGRAD_MAX_ITER) {

		if(gamma >= 0) {
			UMATH_ERROR("minimize_grad", gamma, INVALID_ARGUMENT);
			return vec2(nan());
		}

		vec2 x = guess;
		vec2 grad;
		unsigned int iter = 0;

		do {

			grad = gradient(f, x);
			x += gamma * grad;
			iter++;

		} while(grad.lenght() > MINGRAD_TOLERANCE && iter <= max_iter);

		if(iter > max_iter) {
			UMATH_ERROR("minimize_grad", iter, NO_ALGO_CONVERGENCE);
			return vec2(nan());
		}

		return x;
	}


	/// Find a local maximum of the given multivariate function
	/// using fixed-step gradient descent
	///
	/// @param f The function to minimize
	/// @param guess The initial guess, defaults to {0, 0}
	/// @param gamma The fixed step size, defaults to MINGRAD_GAMMA
	/// @param tolerance The maximum magnitude of the gradient to stop
	/// the algorithm at, defaults to MINGRAD_TOLERANCE.
	/// @param max_iter The maximum number of iterations to perform before
	/// stopping execution of the routine.
	/// @return The coordinates of the local maximum
	inline vec2 maximize_grad(
		multidual<2>(*f)(vec<2, multidual<2>>),
		vec2 guess = {0, 0},
		real gamma = MINGRAD_GAMMA,
		real tolerance = MINGRAD_TOLERANCE,
		unsigned int max_iter = MINGRAD_MAX_ITER) {

		return minimize_grad(f, guess, -gamma, tolerance, max_iter);

	}

}


#endif
