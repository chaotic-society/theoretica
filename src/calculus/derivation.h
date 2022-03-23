#ifndef DERIVATION_UROBORO_H
#define DERIVATION_UROBORO_H

#include "../function.h"
#include "../polynomial/polynomial.h"

namespace uroboro {


	// Derivate a polynomial
	polynomial<real> differentiate_polynomial(polynomial<real> p) {

		polynomial<> Dp;

		for (int i = 1; i < p.size(); ++i) {
			Dp.coeff.push_back(p[i] * i);
		}

		return Dp;
	}


	// Basic derivative approximation
	real approx_derivative(real_function f, real x, real dx = 0) {

		dx = (dx == 0 ? (x / DERIV_PREC) : dx);
		return (f(x + dx) - f(x - dx)) / (2.0 * dx);
	}

}


#endif
