
///
/// @file ortho_polyn.h Orthogonal polynomial bases
///

#ifndef THEORETICA_ORTHO_POLYN_H
#define THEORETICA_ORTHO_POLYN_H

#include "./polynomial.h"
#include "../optimization/roots.h"
#include <functional>


namespace theoretica {


	/// Polynomial sequence recurrence formula type
	/// Used for computing orthogonal polynomial basis elements
	using polyn_recurr_formula
		= std::function<polynomial<real>(polynomial<real>, polynomial<real>, unsigned int)>;


	/// Generate a polynomial basis using a recursion formula
	/// @param P0 First polynomial of the sequence
	/// @param P1 Second polynomial of the sequence
	/// @param f Recursion formula
	/// @param n Degree of the final polynomial
	/// @return The resulting polynomial
	inline polynomial<real> gen_polyn_recurr(
		polynomial<real> P0,
		polynomial<real> P1,
		polyn_recurr_formula f,
		unsigned int n) {

		if(n == 0)
			return P0;

		if(n == 1)
			return P1;

		polynomial<real> Pl;
		for (unsigned int l = 2; l <= n; ++l) {
			
			Pl = f(P0, P1, l);

			P0 = P1;
			P1 = Pl;
		}

		return Pl;
	}


	// Legendre polynomials


	/// Recursion formula for Legendre polynomials
	inline polynomial<real> legendre_polyn_recurr(
		polynomial<real> P0, polynomial<real> P1, unsigned int l) {

		return ((2 * l - 1) * P1 * polynomial<real>({0, 1}) - (l - 1) * P0) / l;
	}


	/// Compute the nth Legendre polynomial
	/// @note The result is not normalized
	inline polynomial<real> legendre_polynomial(unsigned int n) {

		// P0 = 1
		// P1 = x

		return gen_polyn_recurr({1}, {0, 1}, legendre_polyn_recurr, n);
	}


	/// Normalization constant for the nth Legendre polynomial
	inline real legendre_polyn_normalization(unsigned int n) {

		return sqrt((2 * n + 1) / 2.0);
	}


	// Associated Legendre polynomial
	inline std::function<real(real)> assoc_legendre_polynomial(unsigned int l, int m) {

		real K = pow(-1, m);
		polynomial<real> L = legendre_polynomial(l);


		if(m < 0) {
			m = -m;
			K *= pow(-1, m) * fact(l + m) / fact(l - m);
		}

		for (int i = 0; i < m; ++i)
			L = deriv(L);


		if(m % 2 == 0) {

			return [=](real x) {
				return pow(1 - x * x, m / 2) * L(x) / K;
			};

		} else {

			return [=](real x) {
				return sqrt(pow(1 - x * x, m)) * L(x) / K;
			};
		}
	}


	// Associated Legendre polynomial
	inline polynomial<real> assoc_legendre_polynomial_even(unsigned int l, int m) {

		if(m % 2 != 0)
			TH_MATH_ERROR("assoc_legendre_polynomial_even", m, IMPOSSIBLE_OPERATION);

		real K = 1;

		if(m < 0) {
			m = -m;
			K *= fact(l + m) / fact(l - m);
		}

		polynomial<real> L = legendre_polynomial(l);
		polynomial<real> P = ipow(polynomial<real>({1, 0, -1}), m / 2);

		for (int i = 0; i < m; ++i)
			L = deriv(L);

		return L * P / K;
	}


	// Laguerre polynomials


	/// Recursion formula for Laguerre polynomials
	inline polynomial<real> laguerre_polyn_recurr(
		polynomial<real> L0, polynomial<real> L1, unsigned int i) {

		return (polynomial<real>({2 * (real) i - 1, -1}) * L1 - (i - 1) * L0) / i;
	}


	/// Compute the nth Laguerre polynomial
	inline polynomial<real> laguerre_polynomial(unsigned int n) {

		// L0 = 1
		// L1 = 1 - x

		return gen_polyn_recurr({1}, {1, -1}, laguerre_polyn_recurr, n);
	}


	// Generalized Laguerre polynomials


	/// Recursion formula for Generalized Laguerre polynomials
	inline polynomial<real> general_laguerre_polyn_recurr(
		polynomial<real> L0, polynomial<real> L1, real alpha, unsigned int i) {

		return (polynomial<real>({2 * (real) i + alpha - 1, -1}) * L1 - (i + alpha - 1) * L0) / i;
	}


	/// Compute the nth Laguerre polynomial
	inline polynomial<real> general_laguerre_polynomial(real alpha, unsigned int n) {

		// L0 = 1
		// L1 = 1 + alpha - x

		return gen_polyn_recurr(
			{1}, {1 + alpha, -1},
			[alpha](polynomial<real> L0, polynomial<real> L1, unsigned int i) {
				return general_laguerre_polyn_recurr(L0, L1, alpha, i);
			}, n);
	}


	// Hermite polynomials


	/// Recursion formula for Hermite polynomials
	inline polynomial<real> hermite_polyn_recurr(
		polynomial<real> H0, polynomial<real> H1, unsigned int i) {

		return polynomial<real>({0, 2}) * H1 - 2 * (i - 1) * H0;
	}


	/// Compute the nth Hermite polynomial
	/// @note The result is not normalized
	inline polynomial<real> hermite_polynomial(unsigned int n) {

		// H0 = 1
		// H1 = 2x

		return gen_polyn_recurr({1}, {0, 2}, hermite_polyn_recurr, n);
	}


	/// Normalization constant for the nth Hermite polynomial
	inline real hermite_polyn_normalization(unsigned int n) {

		return 1.0 / sqrt((1 << n) * fact(n) * SQRTPI);
	}


	// Chebyshev polynomials


	/// Recursion formula for Chebyshev polynomials
	/// The formula is the same for first and second kind polynomials
	inline polynomial<real> chebyshev_polyn_recurr(
		polynomial<real> T0, polynomial<real> T1, unsigned int i) {

		return polynomial<real>({0, 2}) * T1 - T0;
	}


	/// Compute the nth Chebyshev polynomial of the first kind
	/// @note The result is not normalized
	inline polynomial<real> chebyshev1_polynomial(unsigned int n) {

		// T0 = 1
		// T1 = x

		return gen_polyn_recurr({1}, {0, 1}, chebyshev_polyn_recurr, n);
	}


	/// Compute the nth Chebyshev polynomial of the second kind
	/// @note The result is not normalized
	inline polynomial<real> chebyshev2_polynomial(unsigned int n) {

		// U0 = 1
		// U1 = 2x

		return gen_polyn_recurr({1}, {0, 2}, chebyshev_polyn_recurr, n);
	}


	/// Roots of the n-th Legendre polynomial
	/// @param n The degree of the polynomial
	/// @return A list of the n roots of the Legendre polynomial
	inline std::vector<real> legendre_roots(unsigned int n) {

		if(n == 0)
			return {};

		if(n == 1)
			return {0};

		polynomial<real> P = legendre_polynomial(n);
		std::vector<real> roots;
		roots.reserve(n);

		for (unsigned int i = 1; i <= n; ++i) {
			roots.push_back(
				root_newton(P, (2.0 / (n + 1)) * i - 1.0)
			);
		}

		return roots;
	}


	/// Legendre weights for Gauss-Legendre quadrature of n-th order
	/// @param roots The n roots of the n-th Legendre polynomial
	/// @return The Legendre weights for Gauss-Legendre quadrature
	inline std::vector<real> legendre_weights(const std::vector<real>& roots) {

		const unsigned int n = roots.size();
		const polynomial<real> dP = deriv(legendre_polynomial(n));

		std::vector<real> weights;
		weights.reserve(n);

		for (unsigned int i = 0; i < n; ++i) {
			weights.push_back(
				2.0 / ((1 - square(roots[i])) * square(dP(roots[i])))
			);
		}

		return weights;
	}


	/// Laguerre weights for Gauss-Laguerre quadrature of n-th order
	/// @param roots The n roots of the n-th Laguerre polynomial
	/// @return The Laguerre weights for Gauss-Laguerre quadrature
	inline std::vector<real> laguerre_weights(const std::vector<real>& roots) {

		const unsigned int n = roots.size();
		const polynomial<real> L = laguerre_polynomial(n + 1);

		std::vector<real> weights;
		weights.reserve(n);

		for (unsigned int i = 0; i < n; ++i) {
			weights.push_back(
				roots[i] / square((n + 1) * L(roots[i]))
			);
		}

		return weights;
	}


	/// Hermite weights for Gauss-Hermite quadrature of n-th order
	/// @param roots The n roots of the n-th Hermite polynomial
	/// @return The Hermite weights for Gauss-Hermite quadrature
	inline std::vector<real> hermite_weights(const std::vector<real>& roots) {

		const unsigned int n = roots.size();
		const polynomial<real> H = hermite_polynomial(n - 1);

		std::vector<real> weights;
		weights.reserve(n);

		for (unsigned int i = 0; i < n; ++i) {
			weights.push_back(
				((1 << (n - 1)) * fact(n) * SQRTPI) / square(n * H(roots[i]))
			);
		}

		return weights;
	}

}


#endif
