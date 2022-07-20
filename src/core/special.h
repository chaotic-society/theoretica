
///
/// @file special.h Special functions
///

#ifndef THEORETICA_SPECIAL_H
#define THEORETICA_SPECIAL_H

#include "constants.h"
#include "real_analysis.h"


namespace theoretica {


	/// Compute the real Gamma function
	inline real gamma(real x) {

		real x_fract = fract(x);

		// Identity with the factorial
		if(x_fract <= MACH_EPSILON) {

			if(x >= 1) {
				return fact(int(x) - 1);
			} else {
				UMATH_ERROR("gamma", x, OUT_OF_DOMAIN);
				return nan();
			}
		}

		real x_res = x_fract + 1;
		real mul = 1;

		// Recursion relation for the Gamma function
		// used for domain reduction to [1, 2]
		while(x > 2) {
			x -= 1;
			mul *= x;
		}

		while(x < 1) {
			x += 1;
			mul *= x;
		}

		// Sixth degree interpolating polynomial in [1, 2]
		return mul * (3.0569
			+ x_res * (-4.34693
				+ x_res * (3.25067
					+ x_res * (-1.12613
						+ x_res * 0.165215))));
	}

}

#endif
