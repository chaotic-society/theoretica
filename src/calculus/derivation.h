
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
	template<typename T>
	inline polynomial<T> deriv_polynomial(polynomial<T> p) {

		polynomial<T> Dp;

		for (unsigned int i = 1; i < p.size(); ++i)
			Dp.coeff.push_back(p[i] * i);

		return Dp;
	}


	/// Central derivative approximation using the central method
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param dx The difference in x to use
	inline real deriv_central(real_function f, real x, real dx = DERIV_DX) {
		return (f(x + dx) - f(x - dx)) / (2.0 * dx);
	}


	/// Forward derivative approximation using the forward method
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param dx The difference in x to use
	inline real deriv_forward(real_function f, real x, real dx = DERIV_DX) {
		return (f(x + dx) - f(x)) / dx;
	}


	/// Backward derivative approximation using the backward method
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	/// @param dx The difference in x to use
	inline real deriv_backward(real_function f, real x, real dx = DERIV_DX) {
		return (f(x) - f(x - dx)) / dx;
	}


	/// Use the best available algorithm to approximate
	/// the derivative of a real function
	///
	/// @param f The function to approximate the derivative of
	/// @param x The real value to approximate at
	inline real deriv(real_function f, real x) {
		return deriv_central(f, x);
	}

}


#endif
