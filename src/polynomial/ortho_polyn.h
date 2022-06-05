#ifndef THEORETICA_ORTHO_POLYN_H
#define THEORETICA_ORTHO_POLYN_H

#include "./polynomial.h"


namespace theoretica {


	/// Compute the nth Legendre polynomial
	/// @note The result is not normalized
	inline polynomial<real> legendre_polynomial(unsigned int n) {

		polynomial<real> P0 = {1};
		polynomial<real> P1 = {0, 1};

		if(n == 0)
			return P0;

		if(n == 1)
			return P1;

		polynomial<real> Pl;
		for (int l = 2; l <= n; ++l) {
			
			Pl = ((2 * l - 1) * P1 * polynomial<real>({0, 1}) - (l - 1) * P0) / l;

			P0 = P1;
			P1 = Pl;
		}

		return Pl;
	}


	/// Compute the nth Laguerre polynomial
	inline polynomial<real> laguerre_polynomial(unsigned int n) {

		polynomial<real> L0 = {1};
		polynomial<real> L1 = {1, -1};

		if(n == 0)
			return L0;

		if(n == 1)
			return L1;

		polynomial<real> Li;
		for (int i = 2; i <= n; ++i) {
			
			Li = (polynomial<real>({2 * (real) i - 1, -1}) * L1 - (i - 1) * L0) / i;

			L0 = Li;
			L1 = Li;
		}

		return Li;
	}


	/// Compute the nth Hermite polynomial
	/// @note The result is not normalized
	inline polynomial<real> hermite_polynomial(unsigned int n) {

		polynomial<real> H0 = {1};
		polynomial<real> H1 = {0, 2};

		if(n == 0)
			return H0;

		if(n == 1)
			return H1;

		polynomial<real> Hi;
		for (int i = 2; i <= n; ++i) {
			
			Hi = polynomial<real>({0, 2}) * H1 - 2 * (i - 1) * H0;

			H0 = Hi;
			H1 = Hi;
		}

		return Hi;

	}

}


#endif
