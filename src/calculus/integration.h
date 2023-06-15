
///
/// @file integration.h Integral approximation
///

#ifndef THEORETICA_INTEGRATION_H
#define THEORETICA_INTEGRATION_H

#include "../core/constants.h"
#include "../core/function.h"
#include "../polynomial/polynomial.h"
#include "../polynomial/ortho_polyn.h"
#include "./gauss.h"


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

		res += 0.5 * (f(a) + f(b));

		for (unsigned int i = 1; i < steps; ++i)
			res += f(a + i * dx);

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

		// Sum terms by order of magnitude supposing that
		// f stays at the same order inside the interval,
		// to alleviate truncation errors

		res += f(a) + f(b);

		for (unsigned int i = 2; i < steps; i += 2)
			res += 2 * f(a + i * dx);

		for (unsigned int i = 1; i < steps; i += 2)
			res += 4 * f(a + i * dx);

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
		unsigned int iter = 8) {

		real T[iter][iter];

		T[0][0] = (f(a) + f(b)) * (b - a) / 2.0;

		for (unsigned int j = 1; j < iter; ++j) {
			
			// Composite Trapezoidal Rule
			T[j][0] = integral_trapezoid(f, a, b, 1 << j);

			// Richardson extrapolation
			for (unsigned int k = 1; k <= j; ++k) {
				const uint64_t coeff = 1 << (2 * k);
				T[j][k] = (coeff * T[j][k - 1] - T[j - 1][k - 1]) / (coeff - 1);
			}
		}

		// Return the best approximation
		return T[iter - 1][iter - 1];
	}


	/// Approximate the definite integral of an arbitrary function
	/// using Romberg's method to the given tolerance.
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param tolerance Convergence tolerance for the algorithm
	/// @return An approximation of the integral of f
	template<typename RealFunction>
	inline real integral_romberg_tol(
		RealFunction f,
		real a, real b,
		real tolerance = INTEGRATION_TOL) {

		const unsigned int MAX_ROMBERG_ITER = 16;
		real T[MAX_ROMBERG_ITER][MAX_ROMBERG_ITER];

		T[0][0] = (f(a) + f(b)) * (b - a) / 2.0;

		for (unsigned int j = 1; j < 16; ++j) {
			
			// Composite Trapezoidal Rule
			T[j][0] = integral_trapezoid(f, a, b, 1 << j);

			// Richardson extrapolation
			for (unsigned int k = 1; k <= j; ++k) {
				const uint64_t coeff = 1 << (2 * k);
				T[j][k] = (coeff * T[j][k - 1] - T[j - 1][k - 1]) / (coeff - 1);
			}

			// Stop the algorithm when the desired precision has been reached
			if(abs(T[j][j] - T[j - 1][j - 1]) < tolerance)
				return T[j][j];
		}

		// Return the best approximation
		return T[MAX_ROMBERG_ITER - 1][MAX_ROMBERG_ITER - 1];
	}


	/// Use Gaussian quadrature using the given points and weights.
	///
	/// @param f The function to integrate
	/// @param x The points of evaluation
	/// @param w The weights of the linear combination
	template<typename RealFunction>
	inline real integral_gauss(
		RealFunction f, const std::vector<real>& x, const std::vector<real>& w) {

		if(x.size() != w.size()) {
			TH_MATH_ERROR("integral_gauss", x.size(), INVALID_ARGUMENT);
			return nan();
		}

		real res = 0;

		for (int i = 0; i < x.size(); ++i)
			res += w[i] * f(x[i]);

		return res;
	}


	/// Use Gauss-Legendre quadrature of arbitrary degree to approximate
	/// a definite integral providing the roots of the n degree Legendre polynomial
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param x The roots of the n degree Legendre polynomial
	/// @param w The weights computed for the n-th order quadrature
	/// @return The Gauss-Legendre quadrature of the given function
	template<typename RealFunction>
	inline real integral_legendre(
		RealFunction f, real a, real b, real* x, real* w, unsigned int n) {

		const real mean = (b + a) / 2.0;
		const real halfdiff = (b - a) / 2.0;

		real res = 0;

		for (int i = n - 1; i >= 0; --i)
			res += w[i] * f(halfdiff * x[i] + mean);

		return res * halfdiff;
	}


	/// Use Gauss-Legendre quadrature of arbitrary degree to approximate
	/// a definite integral providing the roots of the n degree Legendre polynomial
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param x The roots of the n degree Legendre polynomial
	/// @param w The weights computed for the n-th order quadrature
	/// @return The Gauss-Legendre quadrature of the given function
	template<typename RealFunction>
	inline real integral_legendre(
		RealFunction f, real a, real b,
		const std::vector<real>& x, const std::vector<real>& w) {

		const real mean = (b + a) / 2.0;
		const real halfdiff = (b - a) / 2.0;

		real res = 0;

		for (int i = x.size() - 1; i >= 0; --i)
			res += w[i] * f(halfdiff * x[i] + mean);

		return res * halfdiff;
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
	inline real integral_legendre(
		RealFunction f, real a, real b, const std::vector<real>& x) {
		
		return integral_legendre(f, a, b, x, legendre_weights(x));
	}


	/// Use Gauss-Legendre quadrature of arbitrary degree to approximate
	/// a definite integral.
	///
	/// @note This function computes the n roots of the n-th degree Legendre
	/// polynomial and the associated weights each time it is called. If multiple
	/// calculations at the same degree are needed, it is more efficient to compute
	/// them only once using legendre_roots.
	///
	/// @note This function uses Newton's method to compute the roots of the
	/// n-th degree Legendre polynomial, so for higher degrees (>> 20) the algorithm
	/// may fail to correctly find all of the zeroes.
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param n The order of the polynomial
	/// @return The Gauss-Legendre quadrature of the given function
	template<typename RealFunction>
	inline real integral_legendre(RealFunction f, real a, real b, unsigned int n = 9) {
		
		switch(n) {
			case 2: return integral_legendre(f, a, b,
				tables::legendre_roots_2, tables::legendre_weights_2, 2); break;
			case 4: return integral_legendre(f, a, b,
				tables::legendre_roots_4, tables::legendre_weights_4, 4); break;
			case 8: return integral_legendre(f, a, b,
				tables::legendre_roots_8, tables::legendre_weights_8, 8); break;
			case 16: return integral_legendre(f, a, b,
				tables::legendre_roots_16, tables::legendre_weights_16, 16); break;
			default: return integral_legendre(f, a, b, legendre_roots(n)); break;
		}
	}


	/// Use Gauss-Laguerre quadrature of arbitrary degree to approximate
	/// an integral over [0, +inf) providing the roots of the n degree Legendre polynomial
	///
	/// @param f The function to integrate
	/// @param x The roots of the n degree Laguerre polynomial
	/// @return The Gauss-Laguerre quadrature of the given function
	template<typename RealFunction>
	inline real integral_laguerre(RealFunction f, const std::vector<real>& x) {
		
		const std::vector<real> weights = laguerre_weights(x);

		real res = 0;

		for (int i = x.size() - 1; i >= 0; --i)
			res += weights[i] * f(x[i]);

		return res;
	}


	/// Use Gauss-Laguerre quadrature of arbitrary degree to approximate
	/// an integral over [a, b] providing the roots of the n degree Legendre polynomial
	///
	/// @param f The function to integrate
	/// @param x The roots of the n degree Laguerre polynomial
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @return The Gauss-Laguerre quadrature of the given function
	template<typename RealFunction>
	inline real integral_laguerre(RealFunction f, real a, real b, const std::vector<real>& x) {
		
		const std::vector<real> weights = laguerre_weights(x);

		const real exp_a = exp(-a);
		const real exp_b = exp(-b);

		real res = 0;

		for (int i = x.size() - 1; i >= 0; --i)
			res += weights[i] * (exp_a * f(x[i] + a) - exp_b * f(x[i] + b));

		return res;
	}


	/// Use Gauss-Hermite quadrature of arbitrary degree to approximate
	/// an integral over (-inf, +inf) providing the roots of the n degree Hermite polynomial
	///
	/// @param f The function to integrate
	/// @param x The roots of the n degree Hermite polynomial
	/// @return The Gauss-Hermite quadrature of the given function
	template<typename RealFunction>
	inline real integral_hermite(RealFunction f, const std::vector<real>& x) {
		
		const std::vector<real> weights = hermite_weights(x);

		real res = 0;

		for (int i = x.size() - 1; i >= 0; --i)
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
		return integral_romberg_tol(f, a, b);
	}


	/// Use the best available algorithm to approximate
	/// the definite integral of a real function to
	/// a given tolerance.
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param tol Tolerance
	/// @return An approximation of the integral of f
	template<typename RealFunction>
	inline real integral(RealFunction f, real a, real b, real tol) {
		return integral_romberg_tol(f, a, b, tol);
	}

}


#endif
