#ifndef DERIVATION_UROBORO_H
#define DERIVATION_UROBORO_H

#include "./function.h"

namespace uroboro {

	constexpr real DERIV_PREC = 100000.0;

	// Basic derivation approximation
	real derivate_approx_base(real_function f, real x, real dx = 0) {

		dx = dx == 0.0 ? (x / DERIV_PREC) : dx;
		return (f(x + dx) - f(x)) / dx;
	}

	// TO-DO Refined approximation of derivative

}


#endif
