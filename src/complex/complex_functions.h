#ifndef UROBORO_COMPLEX_FUNCTIONS
#define UROBORO_COMPLEX_FUNCTIONS

#include "./complex.h"
#include "../common.h"

namespace uroboro {


	// Compute the square of a complex number
	inline complex square(complex z) {
		return complex(
			square(z.Re()) - square(z.Im()),
			2 * z.Re() * z.Im());
	}


	// Compute the cube of a complex number
	inline complex cube(complex z) {
		return complex(
			cube(z.Re()) - 3 * z.Re() * square(z.Im()),
			3 * square(z.Re()) * z.Im() - cube(z.Im()));
	}
	

	// Compute the complex exponential
	inline complex exp(complex z) {
		return complex(cos(z.Im()), sin(z.Im())) * uroboro::exp(z.Re());
	}


	// Return the modulus of a complex number
	inline real abs(complex z) {
		return z.modulus();
	}


	// Computer the complex sine
	inline complex sin(complex z) {

		complex t = z * complex(0, 1);
		return (exp(t) - exp(-t)) / complex(0, 2);
	}


	// Compute the complex cosine
	inline complex cos(complex z) {

		complex t = z * complex(0, 1);
		return (exp(t) + exp(-t)) / 2.0;
	}


	// Compute the complex tangent
	inline complex tan(complex z) {

		complex t = z * complex(0, 2);
		return (exp(t) - 1) / (exp(t) + 1) * complex(0, -1);
	}


	// Compute the complex square root
	inline complex sqrt(complex z) {
		return complex(
			INVSQR2 * sqrt(z.modulus() + z.Re()),
			INVSQR2 * sqrt(z.modulus() - z.Re()));
	}


	// Compute the complex logarithm
	inline complex ln(complex z) {
		return complex(ln(z.modulus()), z.arg());
	}


	// Compute the complex arcsine
	inline complex asin(complex z) {
		return ln(complex(0, 1) * z + sqrt(complex(1, 0) - square(z))) * complex(0, -1);
	}


	// Compute the complex arccosine
	inline complex acos(complex z) {
		return ln(z + sqrt(square(z) - 1)) * complex(0, -1);
	}


	// Compute the complex arctangent
	inline complex atan(complex z) {
		return ln((complex(0, 1) - z) / (complex(0, 1) + z)) * complex(0, -0.5);
	}


}

#endif
