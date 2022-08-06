
///
/// @file derivation.h Derivative approximation
///

#ifndef DERIVATION_THEORETICA_H
#define DERIVATION_THEORETICA_H

#include "../core/function.h"
#include "../polynomial/polynomial.h"


namespace theoretica {


	/// Compute the derivative of a polynomial
	///
	/// @param p The polynomial to differentiate
	/// @return The derivative polynomial
	template<typename T = real>
	inline polynomial<T> deriv_polynomial(const polynomial<T>& p) {

		polynomial<T> Dp;

		for (unsigned int i = 1; i < p.size(); ++i)
			Dp.coeff.push_back(p.get(i) * i);

		return Dp;
	}


	/// Derivative approximation using the central method
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param h The stepsize to use in the finite differences method
	/// @return The approximated value of the derivative
	template<typename RealFunction>
	inline real deriv_central(RealFunction f, real x, real h = DERIV_STEPSIZE) {
		return (f(x + h) - f(x - h)) / (2.0 * h);
	}


	/// Derivative approximation using the forward method
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param h The stepsize to use in the finite differences method
	/// @return The approximated value of the derivative
	template<typename RealFunction>
	inline real deriv_forward(RealFunction f, real x, real h = DERIV_STEPSIZE) {
		return (f(x + h) - f(x)) / h;
	}


	/// Derivative approximation using the backward method
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param h The stepsize to use in the finite differences method
	/// @return The approximated value of the derivative
	template<typename RealFunction>
	inline real deriv_backward(RealFunction f, real x, real h = DERIV_STEPSIZE) {
		return (f(x) - f(x - h)) / h;
	}


	/// Ridder's derivative approximation of second degree
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param h The stepsize to use in the finite differences method
	/// @return The approximated value of the derivative
	template<typename RealFunction>
	inline real deriv_ridders2(RealFunction f, real x, real h = DERIV_STEPSIZE) {
		return (4.0 * deriv_central(f, x, h / 2.0) - deriv_central(f, x, h)) / 3.0;
	}


	/// Ridder's derivative approximation of arbitrary degree
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param degree The degree of the algorithm
	/// @param h The stepsize to use in the finite differences method
	/// @return The approximated value of the derivative
	template<typename RealFunction>
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


	/// Use the best available algorithm to approximate
	/// the derivative of a real function
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param h The stepsize to use in the finite differences method
	/// @return The approximated value of the derivative
	template<typename RealFunction>
	inline real deriv(RealFunction f, real x, real h = DERIV_STEPSIZE) {
		return deriv_ridders2(f, x, h);
	}


	/// Use the best available algorithm to approximate the second
	/// derivative of a real function
	///
	/// @param f The function to approximate the second derivative of
	/// @param x The real value to approximate at
	/// @param h The stepsize to use in the finite differences method
	/// @return The approximated value of the second derivative
	template<typename RealFunction>
	inline real deriv2(RealFunction f, real x, real h = DERIV_STEPSIZE) {
		return (f(x + h) - (2 * f(x)) + f(x - h)) / (h * h);
	}

}


#endif
