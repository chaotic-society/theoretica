
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
#ifdef THEORETICA_HAS_CPP20
	template <
		VectorType Vector,
		VectorType ReturnVector = Vector,
		autodiff::ADScalarField ObjectiveFunction
	>
#else
	template <
		typename Vector = vec<real>,
		typename ReturnVector = Vector,
		typename ObjectiveFunction,
		autodiff::enable_scalar_field<ObjectiveFunction> = true
	>
#endif
	inline ReturnVector multi_minimize_grad(
		ObjectiveFunction f,
		Vector guess,
		real gamma = OPTIMIZATION_MINGRAD_GAMMA,
		real tolerance = OPTIMIZATION_MINGRAD_TOLERANCE,
		unsigned int max_iter = OPTIMIZATION_MINGRAD_ITER
	) {

		ReturnVector x = guess;

		if(gamma >= 0) {
			TH_MATH_ERROR("multi_minimize_grad", gamma, MathError::InvalidArgument);
			return algebra::vec_error(x);
		}

		ReturnVector grad;
		unsigned int iter = 0;

		do {

			grad = gradient(f, x);
			x += gamma * grad;
			iter++;

		} while(grad.norm() > tolerance && iter <= max_iter);

		if(iter > max_iter) {
			TH_MATH_ERROR("multi_minimize_grad", iter, MathError::NoConvergence);
			return algebra::vec_error(x);
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
#ifdef THEORETICA_HAS_CPP20
	template <
		VectorType Vector,
		VectorType ReturnVector = Vector,
		autodiff::ADScalarField ObjectiveFunction
	>
#else
	template <
		typename Vector = vec<real>,
		typename ReturnVector = Vector,
		typename ObjectiveFunction,
		autodiff::enable_scalar_field<ObjectiveFunction> = true
	>
#endif
	inline ReturnVector multi_maximize_grad(
		ObjectiveFunction f,
		Vector guess,
		real gamma = OPTIMIZATION_MINGRAD_GAMMA,
		real tolerance = OPTIMIZATION_MINGRAD_TOLERANCE,
		unsigned int max_iter = OPTIMIZATION_MINGRAD_ITER
	) {
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
#ifdef THEORETICA_HAS_CPP20
	template <
		VectorType Vector,
		VectorType ReturnVector = Vector,
		autodiff::ADScalarField ObjectiveFunction
	>
#else
	template <
		typename Vector = vec<real>,
		typename ReturnVector = Vector,
		typename ObjectiveFunction,
		autodiff::enable_scalar_field<ObjectiveFunction> = true
	>
#endif
	inline ReturnVector multi_minimize_lingrad(
		ObjectiveFunction f,
		Vector guess,
		real tolerance = OPTIMIZATION_MINGRAD_TOLERANCE,
		unsigned int max_iter = OPTIMIZATION_MINGRAD_ITER
	) {

		ReturnVector x = guess;
		ReturnVector grad;
		unsigned int iter = 0;

		constexpr size_t N = ReturnVector::size_argument;

		do {

			grad = autodiff::gradient(f, x);

			// Minimize f(x + gamma * gradient) in [-1, 0]
			// using Golden Section extrema search
			real gamma = minimize_golden(
				[f, x, grad](real gamma){
					return f(multidual<N>::make_argument(x + gamma * grad)).Re();
				}, -1, 0);

			// Fallback value
			if(-gamma <= MACH_EPSILON)
				gamma = OPTIMIZATION_MINGRAD_GAMMA;

			x += gamma * grad;
			iter++;

		} while(grad.norm() > tolerance && iter <= max_iter);

		if(iter > max_iter) {
			TH_MATH_ERROR("multi_minimize_lingrad", iter, MathError::NoConvergence);
			return algebra::vec_error(x);
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
#ifdef THEORETICA_HAS_CPP20
	template <
		VectorType Vector,
		VectorType ReturnVector = Vector,
		autodiff::ADScalarField ObjectiveFunction
	>
#else
	template <
		typename Vector = vec<real>,
		typename ReturnVector = Vector,
		typename ObjectiveFunction,
		autodiff::enable_scalar_field<ObjectiveFunction> = true
	>
#endif
	inline ReturnVector multi_maximize_lingrad(
		ObjectiveFunction f,
		Vector guess,
		real tolerance = OPTIMIZATION_MINGRAD_TOLERANCE,
		unsigned int max_iter = OPTIMIZATION_MINGRAD_ITER) {

		ReturnVector x = guess;
		ReturnVector grad;
		unsigned int iter = 0;

		constexpr size_t N = ReturnVector::size_argument;

		do {

			grad = -gradient(f, x);

			// Maximize f(x + gamma * gradient) in [-1, 0]
			// using Golden Section extrema search
			real gamma = maximize_golden(
				[f, x, grad](real gamma){
					return f(multidual<N>::make_argument(x + gamma * grad)).Re();
				}, -1, 0);

			// Fallback value
			if(-gamma <= MACH_EPSILON)
				gamma = OPTIMIZATION_MINGRAD_GAMMA;

			x += gamma * grad;
			iter++;

		} while(grad.norm() > tolerance && iter <= max_iter);

		if(iter > max_iter) {
			TH_MATH_ERROR("multi_maximize_lingrad", iter, MathError::NoConvergence);
			return algebra::vec_error(x);
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
#ifdef THEORETICA_HAS_CPP20
	template <
		VectorType Vector,
		VectorType ReturnVector = Vector,
		autodiff::ADScalarField ObjectiveFunction
	>
#else
	template <
		typename Vector = vec<real>,
		typename ReturnVector = Vector,
		typename ObjectiveFunction,
		autodiff::enable_scalar_field<ObjectiveFunction> = true
	>
#endif
	inline ReturnVector multi_minimize(
		ObjectiveFunction f, Vector guess,
		real tolerance = OPTIMIZATION_MINGRAD_TOLERANCE
	) {
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
#ifdef THEORETICA_HAS_CPP20
	template <
		VectorType Vector,
		VectorType ReturnVector = Vector,
		autodiff::ADScalarField ObjectiveFunction
	>
#else
	template <
		typename Vector = vec<real>,
		typename ReturnVector = Vector,
		typename ObjectiveFunction,
		autodiff::enable_scalar_field<ObjectiveFunction> = true
	>
#endif
	inline ReturnVector multi_maximize(
		ObjectiveFunction f, Vector guess,
		real tolerance = OPTIMIZATION_MINGRAD_TOLERANCE
	) {
		return multi_maximize_lingrad(f, guess, tolerance);
	}

}


#endif
