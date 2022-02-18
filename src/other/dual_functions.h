#ifndef UROBORO_DUAL_FUNCTIONS_H
#define UROBORO_DUAL_FUNCTIONS_H

#include "./dual.h"
#include "../real_analysis.h"


// Dual numbers can easily be used for
// automatic differentiation, as they behave
// under addition, multiplication and divisions
// as the first derivative.

// Normal operations can be performed and the result
// will have a real part equal to the function evaluated
// for the given argument and a "dual" part equal
// to the first derivative evaluated for the given argument.


namespace uroboro {


	// Return the square of a dual number
	dual square(dual x) {
		return x * x;
	}

	// Return the cube of a dual number
	dual cube(dual x) {
		return x * x * x;
	}

	// Compute the n-th power of a dual number
	dual pow(dual x, int n) {
		real pow_n_1_x = pow(x.Re(), n - 1);
		return dual(pow_n_1_x * x.Re(), pow_n_1_x * n * x.Dual());
	}

	// Compute the square root of a dual number
	dual sqrt(dual x) {
		real sqrt_x = sqrt(x.Re());
		return dual(sqrt_x, 0.5 / sqrt_x * x.Dual());
	}


	dual sin(dual x) {
		return dual(sin(x.Re()), cos(x.Re()) * x.Dual());
	}


	dual cos(dual x) {
		return dual(cos(x.Re()), -sin(x.Re()) * x.Dual());
	}


	dual tan(dual x) {
		return dual(tan(x.Re()), x.Dual() / square(cos(x.Re())));
	}


	dual exp(dual x) {
		return dual(exp(x.Re()), x.Dual() * exp(x.Re()));
	}


	dual ln(dual x) {
		return dual(ln(x.Re()), x.Dual() / x.Re());
	}


	dual abs(dual x) {
		return dual(abs(x.Re()), x.Dual() * sgn(x.Re()));
	}



}


#endif
