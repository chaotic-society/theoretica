
///
/// @file taylor_expansion.h Computation of first and second
/// order Taylor expansions for generic functions using automatic
/// differentiation.
///

#ifndef THEORETICA_TAYLOR_EXPANSION_H
#define THEORETICA_TAYLOR_EXPANSION_H

#include "../polynomial/polynomial.h"
#include "../autodiff/dual.h"
#include "../autodiff/dual2.h"


namespace theoretica {


	/// Computes the first order Taylor expansion of a generic function
	/// around x0, computed using dual numbers.
	inline polynomial<real> taylor_linear_expansion(dual(*f)(dual), real x0 = 0) {

		dual d = f(dual(x0, 1));
		real fx = d.Re();
		real dfx = d.Dual();

		polynomial<real> P = {fx};
		P += polynomial<>({-x0, 1}) * dfx;
		
		return P;
	}


	/// Computes the second order Taylor expansion of a generic function
	/// around x0, computed using dual numbers (of second order).
	inline polynomial<real> taylor_quadratic_expansion(dual2(*f)(dual2), real x0 = 0) {

		dual2 d = f(dual2(x0, 1, 0));
		real fx = d.Re();
		real dfx = d.Dual1();
		real d2fx = d.Dual2();

		polynomial<real> P = {fx};
		P += polynomial<>({-x0, 1}) * dfx;
		P += polynomial<>({square(x0), -2 * x0, 1}) * (d2fx / 2.0);

		return P;
	}

}


#endif
