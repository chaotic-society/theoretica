
///
/// @file complex_analysis.h Functions of complex variable
///

#ifndef THEORETICA_COMPLEX_FUNCTIONS
#define THEORETICA_COMPLEX_FUNCTIONS

#include "./complex.h"
#include "../core/real_analysis.h"


namespace theoretica {


	/// Complex identity
	/// @param z A complex number
	template<typename T>
	inline complex<T> identity(complex<T> z) {
		return z;
	}


	/// Compute the conjugate of a complex number
	/// @param z A complex number
	template<typename T>
	inline complex<T> conjugate(complex<T> z) {
		return complex<T>(z.a, -z.b);
	}


	/// Compute the conjugate of a complex number
	/// @param z A complex number
	template<typename T>
	inline complex<T> inverse(complex<T> z) {
		return conjugate(z) / z.sqr_norm();
	}


	/// Compute the square of a complex number
	/// @param z A complex number
	template<typename T>
	inline complex<T> square(complex<T> z) {
		return complex<T>(
			square(z.Re()) - square(z.Im()),
			2 * z.Re() * z.Im()
		);
	}


	/// Compute the cube of a complex number
	/// @param z A complex number
	template<typename T>
	inline complex<T> cube(complex<T> z) {
		return complex<T>(
			cube(z.Re()) - 3 * z.Re() * square(z.Im()),
			3 * square(z.Re()) * z.Im() - cube(z.Im()));
	}
	

	/// Compute the complex exponential
	/// @param z A complex number
	template<typename T>
	inline complex<T> exp(complex<T> z) {
		return complex<T>(cos(z.Im()), sin(z.Im())) * exp(z.Re());
	}


	/// Return the modulus of a complex number
	/// @param z A complex number
	template<typename T>
	inline real abs(complex<T> z) {
		return z.norm();
	}


	/// Computer the complex sine
	/// @param z A complex number
	template<typename T>
	inline complex<T> sin(complex<T> z) {

		const complex<T> t = z * complex<T>(0, 1);
		return (exp(t) - exp(-t)) * complex<T>(0, -0.5);
	}


	/// Compute the complex cosine
	/// @param z A complex number
	template<typename T>
	inline complex<T> cos(complex<T> z) {

		const complex<T> t = z * complex<T>(0, 1);
		return (exp(t) + exp(-t)) / 2.0;
	}


	/// Compute the complex tangent
	/// @param z A complex number
	template<typename T>
	inline complex<T> tan(complex<T> z) {

		const complex<T> t = exp(z * complex<T>(0, 2));
		return (t - complex<T>(1, 0)) / (t + complex<T>(1, 0)) * complex<T>(0, -1);
	}


	/// Compute the complex square root
	/// @param z A complex number
	template<typename T>
	inline complex<T> sqrt(complex<T> z) {

		if(abs(z.a) < MACH_EPSILON && abs(z.b) < MACH_EPSILON)
			return complex<T>(0);

		return complex<T>(
			INVSQR2 * sqrt((z.norm() + z.Re())),
			INVSQR2 * sqrt((z.norm() - z.Re()))  * sgn(z.b));
	}


	/// Compute the complex logarithm
	/// @param z A complex number
	template<typename T>
	inline complex<T> ln(complex<T> z) {
		return complex<T>(ln(z.norm()), z.arg());
	}


	/// Compute the complex arcsine
	/// @param z A complex number
	template<typename T>
	inline complex<T> asin(complex<T> z) {

		// For real z in [-1, 1], use real asin
		if(abs(z.Im()) < MACH_EPSILON && abs(z.Re()) <= 1.0 + MACH_EPSILON) {
			T real_part = (z.Re() > 1.0) ? 1.0 : ((z.Re() < -1.0) ? -1.0 : z.Re());
			return complex<T>(asin(real_part), 0);
		}

		// General formula with better stability
		return ln(complex<T>(0, 1) * z + sqrt(complex<T>(1, 0) - square(z))) * complex<T>(0, -1);
	}


	/// Compute the complex arccosine
	/// @param z A complex number
	template<typename T>
	inline complex<T> acos(complex<T> z) {

		// For real z in [-1, 1], use real acos
		if(abs(z.Im()) < MACH_EPSILON && abs(z.Re()) <= 1.0 + MACH_EPSILON) {
			const T real_part = (z.Re() > 1.0) ? 1.0 : ((z. Re() < -1.0) ? -1.0 : z.Re());
			return complex<T>(acos(real_part), 0);
		}
		
		// General formula
		return ln(z + sqrt(square(z) - complex<T>(1, 0))) * complex<T>(0, -1);
	}


	/// Compute the complex arctangent
	/// @param z A complex number
	template<typename T>
	inline complex<T> atan(complex<T> z) {
		return ln((complex<T>(0, 1) - z) / (complex<T>(0, 1) + z)) * complex<T>(0, 0.5);
	}


}

#endif
