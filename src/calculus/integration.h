
///
/// @file integration.h Integral approximation
///

#ifndef THEORETICA_INTEGRATION_H
#define THEORETICA_INTEGRATION_H

#include "../core/constants.h"
#include "../core/function.h"
#include "../polynomial/polynomial.h"


namespace theoretica {


	/// Compute the indefinite integral of a polynomial
	///
	/// @param p The polynomial to integrate
	/// @return The indefinite polynomial integral
	template<typename T = real>
	inline polynomial<T> integrate_polynomial(const polynomial<T>& p) {

		polynomial<T> Dp;
		Dp.coeff.push_back(0);

		for (unsigned int i = 0; i < p.size(); ++i)
			Dp.coeff.push_back(p.get(i) / T(i + 1));

		return Dp;
	}


	/// Approximate the definite integral of an arbitrary function
	/// using the midpoint method
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param steps The number of steps
	/// @return An approximation of the integral of f
	template<typename RealFunction>
	inline real integral_midpoint(RealFunction f, real a, real b,
		unsigned int steps = INTEGRATION_STEPS) {

		if(steps == 0) {
			TH_MATH_ERROR("integral_midpoint", steps, DIV_BY_ZERO);
			return nan();
		}
		
		real dx = (b - a) / steps;
		real res = 0;

		for (unsigned int i = 0; i < steps; ++i)
			res += f(a + (i + 0.5) * dx);

		return res * dx;
	}


	/// Approximate the definite integral of an arbitrary function
	/// using the trapezoid method
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param steps The number of steps
	/// @return An approximation of the integral of f
	template<typename RealFunction>
	inline real integral_trapezoid(RealFunction f, real a, real b,
		unsigned int steps = INTEGRATION_STEPS) {

		if(steps == 0) {
			TH_MATH_ERROR("integral_trapezoid", steps, DIV_BY_ZERO);
			return nan();
		}
		
		real dx = (b - a) / steps;
		real res = 0;

		res += 0.5 * f(a);

		for (unsigned int i = 1; i < steps; ++i)
			res += f(a + i * dx);

		res += 0.5 * f(b);

		return res * dx;
	}


	/// Approximate the definite integral of an arbitrary function
	/// using Simpson's method
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param steps The number of steps
	/// @return An approximation of the integral of f
	template<typename RealFunction>
	inline real integral_simpson(RealFunction f, real a, real b,
		unsigned int steps = INTEGRATION_STEPS) {

		if(steps == 0) {
			TH_MATH_ERROR("integral_simpson", steps, DIV_BY_ZERO);
			return nan();
		}
		
		real dx = (b - a) / (real) steps;
		real res = 0;

		res += f(a) + f(b);

		for (unsigned int i = 1; i < steps; ++i) {

			if(i % 2 == 0)
				res += 2.0 * f(a + i * dx);
			else
				res += 4.0 * f(a + i * dx);
			
		}

		return res * dx / 3.0;
	}


	/// Approximate the definite integral of an arbitrary function
	/// using Romberg's method accurate to the given order
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param order The order of accuracy
	/// @return An approximation of the integral of f
	template<typename RealFunction>
	inline real integral_romberg(
		RealFunction f,
		real a, real b,
		unsigned int order = 16) {

		if(order % 2 != 0) {
			TH_MATH_ERROR("integral_romberg", order, IMPOSSIBLE_OPERATION);
			return nan();
		}

		unsigned int iter = order / 2;

		real T[iter][iter];

		T[0][0] = (f(a) + f(b)) * (b - a) / 2.0;

		for (unsigned int j = 1; j < iter; ++j) {
			
			// Composite Trapezoidal Rule
			T[j][0] = integral_trapezoid(f, a, b, 1 << j);

			// Richardson extrapolation
			for (unsigned int k = 1; k <= j; ++k) {
				T[j][k] =
					T[j][k - 1]
					+ (T[j][k - 1] - T[j - 1][k - 1]) / ((1 << (2 * k)) - 1);
			}
		}

		// Return the best approximation
		return T[iter - 1][iter - 1];
	}


	/// Use the best available algorithm to approximate
	/// the definite integral of a real function.
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @return An approximation of the integral of f
	template<typename RealFunction>
	inline real integral(RealFunction f, real a, real b) {
		return integral_simpson(f, a, b);
	}

}


#endif
