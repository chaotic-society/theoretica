
///
/// @file ortho_polyn.h Orthogonal polynomial bases
///

#ifndef THEORETICA_ORTHO_POLYN_H
#define THEORETICA_ORTHO_POLYN_H

#include "./polynomial.h"


namespace theoretica {


	/// Polynomial sequence recurrence formula type
	/// Used for computing orthogonal polynomial basis elements
	using polyn_recurr_formula
		= polynomial<real>(*)(polynomial<real>, polynomial<real>, unsigned int);


	/// Generate a polynomial basis using a recursion formula
	/// @param P0 First polynomial of the sequence
	/// @param P1 Second polynomial of the sequence
	/// @param f Recursion formula
	/// @param 
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

}


#endif
