
///
/// @file derivation.h Derivative approximation
///

#ifndef DERIVATION_THEORETICA_H
#define DERIVATION_THEORETICA_H

#include "../core/function.h"
#include "../polynomial/polynomial.h"


namespace theoretica {


	/// Compute the exact derivative of a polynomial function.
	///
	/// @param p The polynomial to differentiate
	/// @return The derivative polynomial
	template<typename Field = real>
	inline polynomial<Field> deriv(const polynomial<Field>& p) {

		if (p.coeff.size() == 0) {
			TH_MATH_ERROR("deriv", p.coeff.size(), INVALID_ARGUMENT);
			return polynomial<Field>({ static_cast<Field>(nan()) });
		}

		if (p.coeff.size() == 1)
			return polynomial<Field>({static_cast<Field>(0.0)});

		polynomial<Field> Dp;
		Dp.coeff.resize(p.coeff.size() - 1);

		for (unsigned int i = 1; i < p.size(); ++i)
			Dp.coeff[i - 1] = p[i] * i;

		return Dp;
	}


	/// Approximate the first derivative of a real function
	/// using the central method.
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param h The step size to use in the finite differences method
	/// @return The approximated value of the derivative
	template <
		typename RealFunction = std::function<real(real)>,
		enable_real_func<RealFunction> = true
	>
	inline real deriv_central(RealFunction f, real x, real h = CALCULUS_DERIV_STEP) {
		return (f(x + h) - f(x - h)) / (2.0 * h);
	}


	/// Approximate the first derivative of a real function
	/// using the forward method.
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param h The step size to use in the finite differences method
	/// @return The approximated value of the derivative
	template <
		typename RealFunction = std::function<real(real)>,
		enable_real_func<RealFunction> = true
	>
	inline real deriv_forward(RealFunction f, real x, real h = CALCULUS_DERIV_STEP) {
		return (f(x + h) - f(x)) / h;
	}


	/// Approximate the first derivative of a real function
	/// using the backward method.
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param h The step size to use in the finite differences method
	/// @return The approximated value of the derivative
	template <
		typename RealFunction = std::function<real(real)>,
		enable_real_func<RealFunction> = true
	>
	inline real deriv_backward(RealFunction f, real x, real h = CALCULUS_DERIV_STEP) {
		return (f(x) - f(x - h)) / h;
	}


	/// Approximate the first derivative of a real function
	/// using Ridder's method of second degree.
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param h The step size to use in the finite differences method
	/// @return The approximated value of the derivative
	template <
		typename RealFunction = std::function<real(real)>,
		enable_real_func<RealFunction> = true
	>
	inline real deriv_ridders2(RealFunction f, real x, real h = CALCULUS_DERIV_STEP) {
		return (4.0 * deriv_central(f, x, h / 2.0) - deriv_central(f, x, h)) / 3.0;
	}


	/// Approximate the first derivative of a real function 
	/// using Ridder's method of arbitrary degree.
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param degree The degree of the algorithm
	/// @param h The step size to use in the finite differences method
	/// @return The approximated value of the derivative
	template <
		typename RealFunction = std::function<real(real)>,
		enable_real_func<RealFunction> = true
	>
	inline real deriv_ridders(RealFunction f, real x, real h = 0.01, unsigned int degree = 3) {

		real A[degree][degree];

		for (unsigned int i = 0; i < degree; ++i) {
			
			for (unsigned int n = 0; n <= i; ++n) {
				
				unsigned int m = i - n;
				
				if(n == 0) {
					A[n][m] = deriv_central(f, x, h / (1 << m));
				} else {
					real coeff = square(1 << n);
					A[n][m] = (coeff * A[n - 1][m + 1] - A[n - 1][m]) / (coeff - 1);
				}
			}

		}

		return A[degree - 1][0];
	}


	/// Approximate the first derivative of a real function
	/// using the best available algorithm.
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param h The step size to use in the finite differences method
	/// @return The approximated value of the derivative
	template <
		typename RealFunction = std::function<real(real)>,
		enable_real_func<RealFunction> = true
	>
	inline real deriv(RealFunction f, real x, real h = CALCULUS_DERIV_STEP) {
		return deriv_ridders2(f, x, h);
	}


	/// Approximate the second derivative of a real function
	/// using the best available algorithm.
	///
	/// @param f The function to approximate the second derivative of
	/// @param x The real value to approximate at
	/// @param h The step size to use in the finite differences method
	/// @return The approximated value of the second derivative
	template <
		typename RealFunction = std::function<real(real)>,
		enable_real_func<RealFunction> = true
	>
	inline real deriv2(RealFunction f, real x, real h = CALCULUS_DERIV_STEP) {
		return (f(x + h) - (2 * f(x)) + f(x - h)) / (h * h);
	}

}


#endif
