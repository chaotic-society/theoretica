#ifndef UROBORO_APPROXIMATION
#define UROBORO_APPROXIMATION

#include "./constants.h"
#include "./common.h"

namespace uroboro {

	constexpr real APPROXIMATION_TOLERANCE = 0.000001;


	// Approximate e^x using x86 Assembly instructions
	// Works only for positive <x>
	inline real exp_approx(real x) {

		int x_int = abs(int(x - 0.5));
		real x_fract = uroboro::abs(x - x_int);

		// Calculate e^x as e^int(x) * e^fract(x)
		// Where e^fract(x) is calculated as 2^(fract(x) / 2ln2)
		return pow(E, x_int) * square(f2xm1(x_fract / (2 * LN2)) + 1);
	}


	// Approximate powf(real, real) = x^a
	// Using pow(x, int(a)) * exp(fract(a) * ln(x))
	inline real powf_approx(real x, real a) {

		if(a < 0)
			return 1.0 / powf_approx(x, abs(a));

		int a_int = abs(int(a - 0.5));
		real a_fract = uroboro::abs(a - a_int);
		real x_int_pwr = uroboro::pow(x, a_int);

		// Calculate x^fract(a) as e^(x * log2(fract(a) / log2(e)))
		return x_int_pwr * (a_fract >= APPROXIMATION_TOLERANCE ?
			exp_approx(fyl2x(x, a_fract / LOG2E)) : 1);
	}

}


#endif
