
///
/// @file integration.h Integral approximation
///

#ifndef THEORETICA_INTEGRATION_H
#define THEORETICA_INTEGRATION_H

#include "../core/constants.h"
#include "../core/function.h"
#include "../polynomial/polynomial.h"
#include "../polynomial/ortho_polyn.h"


namespace theoretica {


	/// Compute the indefinite integral of a polynomial
	///
	/// @param p The polynomial to integrate
	/// @return The indefinite polynomial integral
	template<typename T = real>
	inline polynomial<T> integrate_polynomial(const polynomial<T>& p) {

		polynomial<T> P;
		P.coeff.push_back(0);

		for (unsigned int i = 0; i < p.size(); ++i)
			P.coeff.push_back(p.get(i) / T(i + 1));

		return P;
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
		
		const real dx = (b - a) / steps;
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
		
		const real dx = (b - a) / steps;
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
		
		const real dx = (b - a) / (real) steps;
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


	/// Use Gauss-Legendre quadrature of arbitrary degree to approximate
	/// a definite integral providing the roots of the n degree Legendre polynomial
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param x The roots of the n degree Legendre polynomial
	/// @return The Gauss-Legendre quadrature of the given function
	template<typename RealFunction>
	inline real integral_gauss_legendre(RealFunction f, real a, real b, const std::vector<real>& x) {
		
		const std::vector<real> weights = legendre_weights(x);
		const real mean = (b + a) / 2.0;
		const real halfdiff = (b - a) / 2.0;

		real res = 0;

		for (unsigned int i = 0; i < x.size(); ++i)
			res += weights[i] * f(halfdiff * x[i] + mean);

		return res * halfdiff;
	}


	/// Use Gauss-Legendre quadrature of arbitrary degree to approximate
	/// a definite integral
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param n The order of the polynomial
	/// @return The Gauss-Legendre quadrature of the given function
	template<typename RealFunction>
	inline real integral_gauss_legendre(RealFunction f, real a, real b, unsigned int n = 9) {
		
		return integral_gauss_legendre(f, a, b, legendre_roots(n));
	}


	/// Use Gauss-Laguerre quadrature of arbitrary degree to approximate
	/// an integral over [0, +inf) providing the roots of the n degree Legendre polynomial
	///
	/// @param f The function to integrate
	/// @param x The roots of the n degree Laguerre polynomial
	/// @return The Gauss-Laguerre quadrature of the given function
	template<typename RealFunction>
	inline real integral_gauss_laguerre(RealFunction f, const std::vector<real>& x) {
		
		const std::vector<real> weights = laguerre_weights(x);

		real res = 0;

		for (unsigned int i = 0; i < x.size(); ++i)
			res += weights[i] * f(x[i]);

		return res;
	}


	/// Use Gauss-Hermite quadrature of arbitrary degree to approximate
	/// an integral over (-inf, +inf) providing the roots of the n degree Hermite polynomial
	///
	/// @param f The function to integrate
	/// @param x The roots of the n degree Hermite polynomial
	/// @return The Gauss-Hermite quadrature of the given function
	template<typename RealFunction>
	inline real integral_gauss_hermite(RealFunction f, const std::vector<real>& x) {
		
		const std::vector<real> weights = hermite_weights(x);

		real res = 0;

		for (unsigned int i = 0; i < x.size(); ++i)
			res += weights[i] * f(x[i]);

		return res;
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
