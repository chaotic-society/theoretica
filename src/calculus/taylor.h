
///
/// @file taylor.h Taylor series expansions
///

#ifndef THEORETICA_TAYLOR_H
#define THEORETICA_TAYLOR_H

#include "../polynomial/polynomial.h"
#include "../autodiff/dual.h"
#include "../autodiff/dual2.h"


namespace theoretica {

	/// @namespace theoretica::taylor Taylor series expansions
	namespace taylor {

		/// Computes the first order Taylor expansion of a generic function
		/// around x0, computed using dual numbers.
		/// Automatic differentiation is used to compute the exact values
		/// of the function and its derivatives at x0.
		///
		/// @param f A function taking and returning a dual number, representing
		/// the function to expand in series.
		/// @param x0 The center of the Taylor expansion
		/// @return The Taylor series expansion of the function to first degree
		template<typename DualFunction>
		inline polynomial<real> linear_expansion(DualFunction f, real x0 = 0) {

			dual d = f(dual(x0, 1));
			real fx = d.Re();
			real dfx = d.Dual();

			polynomial<real> P = {fx};
			P += polynomial<>({-x0, 1}) * dfx;
			
			return P;
		}


		/// Computes the second order Taylor expansion of a generic function
		/// around x0, computed using dual numbers (of second order).
		/// Automatic differentiation is used to compute the exact values
		/// of the function and its derivatives at x0.
		///
		/// @param f A function taking and returning a dual2 number, representing
		/// the function to expand in series.
		/// @param x0 The center of the Taylor expansion
		/// @return The Taylor series expansion of the function to second degree
		template<typename Dual2Function>
		inline polynomial<real> quadratic_expansion(Dual2Function f, real x0 = 0) {

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

}


#endif
