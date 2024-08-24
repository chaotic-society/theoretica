
///
/// @file multi_extrema.h Search of extrema of multivariate functions
///

#ifndef THEORETICA_MULTI_EXTREMA_H
#define THEORETICA_MULTI_EXTREMA_H

#include "../core/constants.h"
#include "../autodiff/autodiff.h"
#include "./extrema.h"

#include <functional>


namespace theoretica {


	/// Find a local minimum of the given multivariate function
	/// using fixed-step gradient descent
	///
	/// @param f The function to multi_minimize
	/// @param guess The initial guess, defaults to the origin
	/// @param gamma The fixed step size, defaults to OPTIMIZATION_MINGRAD_GAMMA
	/// @param tolerance The maximum magnitude of the gradient to stop
	/// the algorithm at, defaults to OPTIMIZATION_MINGRAD_TOLERANCE.
	/// @param max_iter The maximum number of iterations to perform before
	/// stopping execution of the routine.
	/// @return The coordinates of the local minimum, 
	/// NaN if the algorithm did not converge.
	template<unsigned int N>
	inline vec<real, N> multi_minimize_grad(
		multidual<N>(*f)(vec<multidual<N>, N>),
		vec<real, N> guess = vec<real, N>(0),
		real gamma = OPTIMIZATION_MINGRAD_GAMMA,
		real tolerance = OPTIMIZATION_MINGRAD_TOLERANCE,
		unsigned int max_iter = OPTIMIZATION_MINGRAD_ITER) {

		if(gamma >= 0) {
			TH_MATH_ERROR("multi_minimize_grad", gamma, INVALID_ARGUMENT);
			return vec<real, N>(nan());
		}

		vec<real, N> x = guess;
		vec<real, N> grad;
		unsigned int iter = 0;

		do {

			grad = gradient(f, x);
			x += gamma * grad;
			iter++;

		} while(grad.norm() > OPTIMIZATION_MINGRAD_TOLERANCE && iter <= max_iter);

		if(iter > max_iter) {
			TH_MATH_ERROR("multi_minimize_grad", iter, NO_ALGO_CONVERGENCE);
			return vec<real, N>(nan());
		}

		return x;
	}


	/// Find a local maximum of the given multivariate function
	/// using fixed-step gradient descent
	///
	/// @param f The function to multi_maximize
	/// @param guess The initial guess, defaults to the origin
	/// @param gamma The fixed step size, defaults to OPTIMIZATION_MINGRAD_GAMMA
	/// @param tolerance The maximum magnitude of the gradient to stop
	/// the algorithm at, defaults to OPTIMIZATION_MINGRAD_TOLERANCE.
	/// @param max_iter The maximum number of iterations to perform before
	/// stopping execution of the routine.
	/// @return The coordinates of the local maximum, 
	/// NaN if the algorithm did not converge.
	template<unsigned int N>
	inline vec<real, N> multi_maximize_grad(
		multidual<N>(*f)(vec<multidual<N>, N>),
		vec<real, N> guess = vec<real, N>(0),
		real gamma = OPTIMIZATION_MINGRAD_GAMMA,
		real tolerance = OPTIMIZATION_MINGRAD_TOLERANCE,
		unsigned int max_iter = OPTIMIZATION_MINGRAD_ITER) {

		return multi_minimize_grad(f, guess, -gamma, tolerance, max_iter);
	}


	/// Find a local minimum of the given multivariate function
	/// using gradient descent with linear search
	///
	/// @param f The function to multi_minimize
	/// @param guess The initial guess, defaults to the origin
	/// @param tolerance The maximum magnitude of the gradient to stop
	/// the algorithm at, defaults to OPTIMIZATION_MINGRAD_TOLERANCE.
	/// @param max_iter The maximum number of iterations to perform before
	/// stopping execution of the routine.
	/// @return The coordinates of the local minimum, 
	/// NaN if the algorithm did not converge.
	template<unsigned int N>
	inline vec<real, N> multi_minimize_lingrad(
		multidual<N>(*f)(vec<multidual<N>, N>),
		vec<real, N> guess = vec<real, N>(0),
		real tolerance = OPTIMIZATION_MINGRAD_TOLERANCE,
		unsigned int max_iter = OPTIMIZATION_MINGRAD_ITER) {

		vec<real, N> x = guess;
		vec<real, N> grad;
		unsigned int iter = 0;

		do {

			grad = autodiff::gradient(f, x);

			// Minimize f(x + gamma * gradient) in [-1, 0]
			// using Golden Section extrema search
			real gamma = minimize_goldensection(
				[f, x, grad](real gamma){
					return f(
						multidual<N>::make_argument(x)
						+ multidual<N>::make_argument(grad)
						* multidual<N>(gamma)).Re();
				},
				-1, 0);

			// Fallback value
			if(-gamma <= MACH_EPSILON)
				gamma = OPTIMIZATION_MINGRAD_GAMMA;

			x += gamma * grad;
			iter++;

		} while(grad.norm() > OPTIMIZATION_MINGRAD_TOLERANCE && iter <= max_iter);

		if(iter > max_iter) {
			TH_MATH_ERROR("multi_minimize_lingrad", iter, NO_ALGO_CONVERGENCE);
			return vec<real, N>(nan());
		}

		return x;
	}


	/// Find a local maximum of the given multivariate function
	/// using gradient descent with linear search
	///
	/// @param f The function to multi_maximize
	/// @param guess The initial guess, defaults to the origin
	/// @param tolerance The minimum magnitude of the gradient to stop
	/// the algorithm at, defaults to OPTIMIZATION_MINGRAD_TOLERANCE.
	/// @param max_iter The maximum number of iterations to perform before
	/// stopping execution of the routine.
	/// @return The coordinates of the local maximum, 
	/// NaN if the algorithm did not converge.
	template<unsigned int N>
	inline vec<real, N> multi_maximize_lingrad(
		multidual<N>(*f)(vec<multidual<N>, N>),
		vec<real, N> guess = vec<real, N>(0),
		real tolerance = OPTIMIZATION_MINGRAD_TOLERANCE,
		unsigned int max_iter = OPTIMIZATION_MINGRAD_ITER) {

		vec<real, N> x = guess;
		vec<real, N> grad;
		unsigned int iter = 0;

		do {

			grad = -gradient(f, x);

			// Maximize f(x + gamma * gradient) in [-1, 0]
			// using Golden Section extrema search
			real gamma = maximize_goldensection(
				[f, x, grad](real gamma){
					return f(
						multidual<N>::make_argument(x)
						+ multidual<N>::make_argument(grad) * gamma).Re();
				},
				-1, 0);

			// Fallback value
			if(-gamma <= MACH_EPSILON)
				gamma = OPTIMIZATION_MINGRAD_GAMMA;

			x += gamma * grad;
			iter++;

		} while(grad.norm() > OPTIMIZATION_MINGRAD_TOLERANCE && iter <= max_iter);

		if(iter > max_iter) {
			TH_MATH_ERROR("multi_maximize_lingrad", iter, NO_ALGO_CONVERGENCE);
			return vec<real, N>(nan());
		}

		return x;
	}


	/// Use the best available algorithm to find a local
	/// minimum of the given multivariate function
	///
	/// @param f The function to multi_minimize
	/// @param guess The initial guess
	/// @param tolerance The minimum magnitude of the gradient to stop
	/// the algorithm at, defaults to OPTIMIZATION_MINGRAD_TOLERANCE.
	/// @return The coordinates of the local minimum, 
	/// NaN if the algorithm did not converge.
	template<unsigned int N>
	inline vec<real, N> multi_minimize(
		multidual<N>(*f)(vec<multidual<N>, N>),
		vec<real, N> guess = vec<real, N>(0), real tolerance = OPTIMIZATION_MINGRAD_TOLERANCE) {

		return multi_minimize_lingrad(f, guess, tolerance);
	}


	/// Use the best available algorithm to find a local
	/// maximum of the given multivariate function
	///
	/// @param f The function to multi_maximize
	/// @param guess The initial guess
	/// @param tolerance The minimum magnitude of the gradient to stop
	/// the algorithm at, defaults to OPTIMIZATION_MINGRAD_TOLERANCE.
	/// @return The coordinates of the local maximum, 
	/// NaN if the algorithm did not converge.
	template<unsigned int N>
	inline vec<real, N> multi_maximize(
		multidual<N>(*f)(vec<multidual<N>, N>),
		vec<real, N> guess = vec<real, N>(0), real tolerance = OPTIMIZATION_MINGRAD_TOLERANCE) {

		return multi_maximize_lingrad<N>(f, guess, tolerance);
	}

}


#endif
