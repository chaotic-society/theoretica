#ifndef UROBORO_APPROXIMATION
#define UROBORO_APPROXIMATION

#include "./constants.h"
#include "./common.h"

namespace uroboro {

	constexpr real APPROXIMATION_TOLERANCE = 0.000001;


	// Approximate e^x using power series
	inline real exp_approx(real x) {

		int x_int = int(x - 0.5);
		real x_fract = uroboro::abs(x - x_int);

		return pow(E, x_int) * square(f2xm1(x_fract / (2 * LN2)) + 1);
	}


	// Approximate powf(real, real) = x^a
	// Using pow(x, int(a)) * exp(fract(a) * ln(x))
	// Works only if <a> is positive
	inline real powf_approx(real x, real a) {

		int a_int = abs(int(a - 0.5));
		real a_fract = uroboro::abs(a - a_int);
		real x_int_pwr = uroboro::pow(x, a_int);

		if(a < 0)
			x_int_pwr = 1.0 / x_int_pwr;


		if(a_fract >= APPROXIMATION_TOLERANCE)
			return x_int_pwr * exp_approx(fyl2x(x, a_fract / LOG2E));
		else
			return x_int_pwr;
	}

}


#endif
