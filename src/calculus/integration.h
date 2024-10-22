
///
/// @file integration.h Integral approximation
///

#ifndef THEORETICA_INTEGRATION_H
#define THEORETICA_INTEGRATION_H

#include "../core/constants.h"
#include "../core/function.h"
#include "../polynomial/polynomial.h"
#include "../polynomial/orthogonal.h"
#include "./gauss.h"


namespace theoretica {


	/// Compute the indefinite integral of a polynomial
	///
	/// @param p The polynomial to integrate
	/// @return The indefinite polynomial integral
	template<typename T = real>
	inline polynomial<T> integral(const polynomial<T>& p) {

		polynomial<T> P;
		P.coeff.resize(p.size() + 1);
		P[0] = 0;

		for (unsigned int i = 0; i < p.size(); ++i)
			P[i + 1] = p[i] / T(i + 1);

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
		unsigned int steps = CALCULUS_INTEGRAL_STEPS) {

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
		unsigned int steps = CALCULUS_INTEGRAL_STEPS) {

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
		unsigned int steps = CALCULUS_INTEGRAL_STEPS) {

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
				const uint64_t coeff = uint64_t(1) << (2 * k);
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
		real tolerance = CALCULUS_INTEGRAL_TOL) {

		const unsigned int MAX_ROMBERG_ITER = 16;
		real T[MAX_ROMBERG_ITER][MAX_ROMBERG_ITER];

		T[0][0] = (f(a) + f(b)) * (b - a) / 2.0;

		for (unsigned int j = 1; j < 16; ++j) {
			
			// Composite Trapezoidal Rule
			T[j][0] = integral_trapezoid(f, a, b, 1 << j);

			// Richardson extrapolation
			for (unsigned int k = 1; k <= j; ++k) {
				const uint64_t coeff = uint64_t(1) << (2 * k);
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


	/// Use Gaussian quadrature using the given points and weights.
	///
	/// @param f The function to integrate
	/// @param x The points of evaluation
	/// @param w The weights of the linear combination
	/// @param n The number of points used
	template<typename RealFunction>
	inline real integral_gauss(
		RealFunction f, real* x, real* w, unsigned int n) {

		real res = 0;

		for (unsigned int i = 0; i < n; ++i)
			res += w[i] * f(x[i]);

		return res;
	}


	/// Use Gaussian quadrature using the given points and weights.
	///
	/// @param f The function to integrate
	/// @param x The points of evaluation
	/// @param w The weights of the linear combination
	/// @param n The number of points used
	/// @param Winv The inverse of the weight function
	template<typename RealFunction>
	inline real integral_gauss(
		RealFunction f, real* x, real* w, unsigned int n,
		real_function Winv) {

		real res = 0;

		for (unsigned int i = 0; i < n; ++i)
			res += w[i] * f(x[i]) * Winv(x[i]);

		return res;
	}


	/// Use Gauss-Legendre quadrature of arbitrary degree to approximate
	/// a definite integral providing the roots of the n degree Legendre polynomial.
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
	/// a definite integral providing the roots of the n degree Legendre polynomial.
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
	/// a definite integral providing the roots of the n degree Legendre polynomial.
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


	/// Use Gauss-Legendre quadrature of degree 2, 4, 8 or 16,
	/// using pre-computed values, to approximate
	/// an integral over [a, b].
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param n The order of the polynomial (available values are
	/// 2, 4, 8, 16 or 32).
	/// @return The Gauss-Legendre quadrature of the given function
	template<typename RealFunction>
	inline real integral_legendre(RealFunction f, real a, real b, unsigned int n = 16) {
		
		switch(n) {
			case 2: return integral_legendre(f, a, b,
				tables::legendre_roots_2, tables::legendre_weights_2, 2); break;
			case 4: return integral_legendre(f, a, b,
				tables::legendre_roots_4, tables::legendre_weights_4, 4); break;
			case 8: return integral_legendre(f, a, b,
				tables::legendre_roots_8, tables::legendre_weights_8, 8); break;
			case 16: return integral_legendre(f, a, b,
				tables::legendre_roots_16, tables::legendre_weights_16, 16); break;
			// case 32: return integral_legendre(f, a, b,
			// 	tables::legendre_roots_32, tables::legendre_weights_32, 32); break;
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

		return integral_gauss(f, x, laguerre_weights(x));
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
	inline real integral_laguerre(
		RealFunction f, real a, real b, const std::vector<real>& x) {
		
		const std::vector<real> weights = laguerre_weights(x);

		const real exp_a = exp(-a);
		const real exp_b = exp(-b);

		real res = 0;

		for (int i = x.size() - 1; i >= 0; --i)
			res += weights[i] * (exp_a * f(x[i] + a) - exp_b * f(x[i] + b));

		return res;
	}


	/// Use Gauss-Laguerre quadrature of degree 2, 4, 8 or 16,
	/// using pre-computed values, to approximate
	/// an integral over [0, +inf).
	///
	/// @param f The function to integrate
	/// @param n The order of the polynomial (available values are
	/// 2, 4, 8 or 16).
	/// @return The Gauss-Legendre quadrature of the given function
	template<typename RealFunction>
	inline real integral_laguerre(RealFunction f, unsigned int n = 16) {
		
		switch(n) {
			case 2: return integral_gauss(f,
				tables::laguerre_roots_2, tables::laguerre_weights_2, 2); break;
			case 4: return integral_gauss(f,
				tables::laguerre_roots_4, tables::laguerre_weights_4, 4); break;
			case 8: return integral_gauss(f,
				tables::laguerre_roots_8, tables::laguerre_weights_8, 8); break;
			case 16: return integral_gauss(f,
				tables::laguerre_roots_16, tables::laguerre_weights_16, 16); break;
			// case 32: return integral_gauss(f,
			// 	tables::laguerre_roots_32, tables::laguerre_weights_32, 32); break;
			default: {
				TH_MATH_ERROR("integral_laguerre", n, INVALID_ARGUMENT);
				return nan(); break;
			}
		}
	}


	/// Use Gauss-Hermite quadrature of arbitrary degree to approximate
	/// an integral over (-inf, +inf) providing the roots of the n degree Hermite polynomial
	///
	/// @param f The function to integrate
	/// @param x The roots of the n degree Hermite polynomial
	/// @return The Gauss-Hermite quadrature of the given function
	template<typename RealFunction>
	inline real integral_hermite(RealFunction f, const std::vector<real>& x) {
		
		return integral_gauss(f, x, hermite_weights(x));
	}


	/// Use Gauss-Hermite quadrature of degree 2, 4, 8 or 16,
	/// using pre-computed values, to approximate
	/// an integral over (-inf, +inf).
	///
	/// @param f The function to integrate
	/// @param n The order of the polynomial (available values are
	/// 2, 4, 8 or 16).
	/// @return The Gauss-Hermite quadrature of the given function
	template<typename RealFunction>
	inline real integral_hermite(RealFunction f, unsigned int n = 16) {
		
		switch(n) {
			case 2: return integral_gauss(f,
				tables::hermite_roots_2, tables::hermite_weights_2, 2); break;
			case 4: return integral_gauss(f,
				tables::hermite_roots_4, tables::hermite_weights_4, 4); break;
			case 8: return integral_gauss(f,
				tables::hermite_roots_8, tables::hermite_weights_8, 8); break;
			case 16: return integral_gauss(f,
				tables::hermite_roots_16, tables::hermite_weights_16, 16); break;
			default: {
				TH_MATH_ERROR("integral_hermite", n, INVALID_ARGUMENT);
				return nan(); break;
			}
		}
	}


	/// Integrate a function from a point up to infinity
	/// by integrating it by steps, stopping execution when
	/// the variation of the integral is small enough or
	/// the number of steps reaches a maximum value.
	inline real integral_inf_riemann(
		real_function f, real a, real step_sz = 1,
		real tol = CALCULUS_INTEGRAL_TOL, unsigned int max_iter = 100) {

		// Current lower extreme of the interval
		real x_n = a + step_sz;

		// Total integral sum
		real sum = integral_romberg_tol(f, a, x_n, tol);
		
		// Variation between steps
		real delta = inf();

		// Number of steps performed
		unsigned int i = 0;

		while(abs(delta) > tol && i <= max_iter) {

			delta = integral_romberg_tol(f, x_n, x_n + step_sz, tol);
			sum += delta;
			x_n += step_sz;
			i++;
		}

		if(i >= max_iter) {
			TH_MATH_ERROR("integral_inf_riemann", i, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return sum;
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
