
///
/// @file real_analysis.h Real functions
///

#ifndef UROBORO_COMMON_H
#define UROBORO_COMMON_H

#include "./constants.h"
#include "./error.h"


namespace uroboro {


	/// Compute the square of a real number
	/// @param[in] x A real number
	/// @param[out] The square of <x>
	inline real square(real x) {
		return x * x;
	}


	/// Compute the cube of a real number
	/// @param[in] x A real number
	/// @param[out] The cube of <x>
	inline real cube(real x) {
		return x * x * x;
	}


	/// Compute the square root of a real number
	/// @param[in] x A real number
	/// @param[out] The square root of x 
	inline real sqrt(real x) {

		if(x < 0) {
			UMATH_ERROR("sqrt", x, OUT_OF_DOMAIN);
			return nan();
		}

#ifdef UROBORO_X86

	#ifdef MSVC_ASM
		__asm {
			fld x
			fsqrt
		}
	#else
		asm("fsqrt" : "+t"(x));
	#endif

		return x;
#else

		if(x < 1) {

			if(x == 0)
				return 0;

			// Approximate sqrt(x) between 0 and 1
			// The root of the inverse is the inverse of the root
			// !!! Possible precision problems with smaller numbers
			return 1.0 / sqrt(1.0 / x);
		}

		// Approximate sqrt(x) using Newton-Raphson
		real y = x;
		int i = 0;

		while((square(y) - x) > ROOT_APPROX_TOL && i < MAX_NEWTON_ITER) {
			y = (y + x / y) / 2.0;
			i++;
		}

		return y;
#endif

	}


	/// Compute the cubic root of x
	/// @param[in] x A real number
	/// @param[out] The cubic root of <x>
	inline real cbrt(real x) {

		if(x < 1) {

			if(x == 0)
				return 0;

			// cbrt(x) is odd
			if(x < 0)
				return -cbrt(-x);

			// Approximate cbrt between 0 and 1
			// The root of the inverse is the inverse of the root
			// !!! Possible precision problems with smaller numbers
			return 1.0 / cbrt(1.0 / x);
		}

		// Approximate cbrt(x) using Newton-Raphson

		real y = x;
		int i = 0;

		while((cube(y) - x) > ROOT_APPROX_TOL && i < MAX_NEWTON_ITER) {
			y = (y * 2.0 + x / (y * y)) / 3.0;
			i++;
		}

		return y;
	}


	/// Compute the absolute value of a real number
	/// @param[in] x A real number
	/// @param[out] The absolute value of <x>
	inline real abs(real x) {

#ifdef UROBORO_X86

	#ifdef MSVC_ASM
		__asm {
			fld x
			fabs
		}
	#else
		asm("fabs" : "+t"(x));
	#endif

		return x;
#else
		return x >= 0 ? x : -x;
#endif

	}


	/// Return the sign of x (1 if positive, -1 if negative, 0 if null)
	/// @param[in] x A real number
	/// @param[out] The sign of <x>
	inline int sgn(real x) {
		return (x > 0) ? 1 : (x < 0 ? -1 : 0);
	}


	/// Compute the floor of x
	/// Computes the maximum integer number that is smaller than x
	/// @param[in] x A real number
	/// @param[out] The floor of <x>
	inline int floor(real x) {

		if(x < 0 && x > -1)
			return -1;

		// Compute the biggest integer number
		// that is smaller than x
		return x - (int(x) % 1);
	}


	/// Compute the fractional part of a real number
	/// @param[in] x A real number
	/// @param[out] The fractional part of <x>
	inline real fract(real x) {
		return abs(x - floor(x));
	}


	/// Return the greatest number between two real numbers
	/// @param[in] x A real number
	/// @param[in] y A real number
	/// @param[out] The greatest number between <x> and <y>
	inline real max(real x, real y) {
		
		#ifdef UROBORO_BRANCHLESS
			return (x + y + abs(x - y)) / 2.0;
		#else
			return x > y ? x : y;
		#endif
	}


	/// Return the smallest number between two real numbers
	/// @param[in] x A real number
	/// @param[in] y A real number
	/// @param[out] The smallest number between <x> and <y>
	inline real min(real x, real y) {

		#ifdef UROBORO_BRANCHLESS
			return (x + y - abs(x - y)) / 2.0;
		#else
			return x > y ? y : x;
		#endif
	}


	/// Clamp x between a and b
	/// @param[in] x A real number
	/// @param[in] a A real number
	/// @param[in] b A real number
	/// @param[out] Returns <x> if <x> is between <a> and <b>, \
	/// <a> if <x> is less than <a>, <b> if <x> is bigger than <b>
	inline real clamp(real x, real a, real b) {

#ifdef UROBORO_FORCE_BRANCHLESS

		// The branchless implementation might be slower or equal
		// on most compilers
		return min(max(x, a), b);
#else
		return x > b ? b : (x < a ? a : x);
#endif
	}


	// x86 instruction wrappers

#ifdef UROBORO_X86

	/// Compute y * log2(x) using x86 Assembly instructions
	inline real fyl2x(real x, real y) {

		real res;

		#ifdef MSVC_ASM

		// TO-DO

		#else
		asm ("fyl2x" : "=t"(res) : "0"(x), "u"(y) : "st(1)");
		#endif

		return res;
	}


	/// Compute 2^x - 1 using x86 Assembly instructions
	/// Works only between -1 and +1
	/// May become particularly incorrect near boundaries
	inline real f2xm1(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			f2xm1
		}
		#else
		asm("f2xm1" : "+t"(x));
		#endif

		return x;
	}

#endif

	/// Compute the binary logarithm of a real number
	/// @param[in] x A real number bigger than 0
	/// @param[out] The binary logarithm of <x>
	inline real log2(real x) {

		if(x <= 0) {

			if(x == 0) {
				UMATH_ERROR("log2", x, OUT_OF_RANGE);
				return -inf();
			}

			UMATH_ERROR("log2", x, OUT_OF_DOMAIN);
			return nan();
		}

#ifdef UROBORO_X86
		return fyl2x(x, 1.0);
#else
		return nan();
#endif
	}


	/// Compute the base-10 logarithm of x
	/// @param[in] x A real number bigger than 0
	/// @param[out] The base-10 logarithm of <x>
	inline real log10(real x) {

		if(x <= 0) {

			if(x == 0) {
				UMATH_ERROR("log10", x, OUT_OF_RANGE);
				return -inf();
			}

			UMATH_ERROR("log10", x, OUT_OF_DOMAIN);
			return nan();
		}

#ifdef UROBORO_X86
		return fyl2x(x, 1.0 / LOG210);
#else
		return nan();
#endif
	}


	/// Compute the natural logarithm of x
	/// @param[in] x A real number bigger than 0
	/// @param[out] The natural logarithm of <x>
	inline real ln(real x) {

		if(x <= 0) {

			if(x == 0) {
				UMATH_ERROR("ln", x, OUT_OF_RANGE);
				return -inf();
			}

			UMATH_ERROR("ln", x, OUT_OF_DOMAIN);
			return nan();
		}

#ifdef UROBORO_X86
		return fyl2x(x, 1.0 / LOG2E);
#else
		return nan();
#endif
	}


	/// Compute the <n>-th power of <x> (where n is natural)
	/// @param[in] x A real number
	/// @param[in] n The integer exponent
	/// @param[out] <x> to the power <n>
	template<typename T>
	inline T pow(T x, int n) {

		T res;
		T x2 = x * x;
		int i = 1;

		if(n > 0) {

			res = x;

			for(; i < (n / 2); i += 2)
				res *= x2;

			for(; i < n; ++i)
				res *= x;

		} else if(n < 0) {

			res = 1 / x;

			for(; i < (n / 2); i += 2)
				res /= x2;

			for(; i < -n; ++i)
				res /= x;
		} else {
			return 1.0;
		}

		return res;
	}


	/// Compute the n-th root of x
	/// @param[in] x A real number
	/// @param[in] n The root number
	/// @param[out] The <n>-th real root of <x>
	inline real root(real x, int n) {

		if(n % 2 == 0 && x < 0 || n == 0) {
			UMATH_ERROR("root", n, OUT_OF_DOMAIN);
			return nan();
		}

		if(n < 0)
			return 1.0 / root(x, -n);

		if(x < 1) {

			if(x == 0)
				return 0;

			// Approximate root between 0 and 1
			// The root of the inverse is the inverse of the root
			// !!! Possible precision problems with smaller numbers
			return 1.0 / root(1.0 / x, n);
		}

		// Approximate n-th root using Newton-Raphson

		real y = x;
		int i = 0;

		while((pow(y, n) - x) > ROOT_APPROX_TOL && i < MAX_NEWTON_ITER) {
			y = (y * (n - 1) + x / pow(y, n - 1)) / (real) n;
			i++;
		}

		return y;
	}


	/// Compute the factorial of n
	inline long long fact(unsigned int n) {

		long long res = 1;
		for (int i = n; i > 1; --i)
			res *= i;

		return res;
	}


#ifdef UROBORO_X86

	/// Approximate e^x using x86 Assembly instructions
	/// Works only for positive <x> in [0, 1]
	inline real exp_x86_norm(real x) {

		// e^x is calculated as 2^(x / ln2)
		return square(f2xm1(x / (2 * LN2)) + 1);
	}

#endif


	/// Compute e^x
	inline real exp(real x) {

	// Domain reduction to [0, +inf]
	if(x < 0)
		return 1.0 / exp(abs(x));

	real fract_x = fract(x);

#ifdef UROBORO_X86

		return pow(E, floor(x)) * exp_x86_norm(fract_x);

#else

	// Taylor series expansion
	// Compute e^floor(x) * e^fract(x)
	
	real res = 1;
	real s_n = 1;

	for (int i = 1; i < TAYLOR_ORDER; ++i) {

		// Recurrence formula to improve
		// numerical stability and performance
		s_n *= fract_x / static_cast<real>(i);
		res += s_n;
	}

	return pow(E, floor(x)) * res;

#endif
	}


	/// Approximate powf(real, real) = x^a
	/// @param[in] x A real number
	/// @param[out] a A real exponent
	inline real powf(real x, real a) {

		if(a < 0)
			return 1.0 / exp(abs(a) * ln(x));

		// x^a = e^(a * ln(x))
		return exp(a * ln(abs(x)));
	}


	/// Compute the sine of a real number
	/// @param[in] x An angle in **radians**
	/// @param[out] The sine of <x>
	inline real sin(real x) {

#ifdef UROBORO_X86

	#ifdef MSVC_ASM
		__asm {
			fld x
			fsin
		}
	#else
		asm("fsin" : "+t"(x));
	#endif

		return x;

#else

		// Clamp x between -2PI and 2PI
		while(x >= 2 * PI)
			x -= 2 * PI;
		
		while(x <= -2 * PI)
			x += 2 * PI;

		real res = 0;

		// Taylor series expansion
		// sin(x) = sum( (-1)^i * x^(2i+1) / (2i+1)! )

		for (int i = 0; i < TAYLOR_ORDER; ++i) {
			res += (i % 2 == 0 ? 1 : -1)
				* pow(x, 2 * i + 1)
				/ static_cast<real>(fact(2 * i + 1));
		}

		return res;
#endif
	}


	/// Compute the cosine of a real number
	/// @param[in] x An angle in **radians**
	/// @param[out] The cosine of <x>
	inline real cos(real x) {

#ifdef UROBORO_X86

		#ifdef MSVC_ASM
		__asm {
			fld x
			fcos
		}
		#else
		asm("fcos" : "+t"(x));
		#endif

		return x;

#else

		// cos(x) is even (cos(x) = cos(-x))
		if(x < 0)
			x = -x;

		// Clamp x between 0 and 2PI
		while(x >= 2 * PI)
			x -= 2 * PI;

		real res = 0;

		// Taylor series expansion
		// sin(x) = sum( (-1)^i * x^(2i) / (2i)! )

		for (int i = 0; i < TAYLOR_ORDER; ++i) {
			res += (i % 2 == 0 ? 1 : -1)
				* pow(x, 2 * i)
				/ static_cast<real>(fact(2 * i));
		}

		return res;
#endif
	}


	/// Compute the tangent of x
	inline real tan(real x) {

		real s, c;

#ifdef UROBORO_X86

		#ifdef MSVC_ASM

		// TO-DO

		#else
		asm ("fsincos" : "=t"(c), "=u"(s) : "0"(x));
		#endif

#else
		s = sin(x);
		c = cos(x);
#endif

		if(c == 0) {
			UMATH_ERROR("tan", c, DIV_BY_ZERO);
		}

		return s / c;
	}


	/// Compute the cotangent of x
	inline real cot(real x) {

		real s, c;

#ifdef UROBORO_X86

		#ifdef MSVC_ASM

		// TO-DO

		#else
		asm ("fsincos" : "=t"(c), "=u"(s) : "0"(x));
		#endif

#else
		s = sin(x);
		c = cos(x);
#endif

		if(s == 0) {
			UMATH_ERROR("cot", s, DIV_BY_ZERO);
		}

		return c / s;
	}


	/// Compute the arctangent
	inline real atan(real x) {

		// Domain reduction to [-1, 1]
		if(abs(x) > 1.0)
			return (PI2 - atan(1.0 / abs(x))) * sgn(x);

		real x2 = x * x;

		// Interpolating Chebyshev polynomial
		// of order 9
		return x * (0.999965
			+ x2 * (-0.331545
				+ x2 * (0.184464
					+ x2 * (-0.090752
						+ 0.023286 * x2))));
	}


	/// Compute the arcsine
	inline real asin(real x) {

		if(abs(x) > 1) {
			UMATH_ERROR("asin", x, OUT_OF_DOMAIN);
			return nan();
		}

		return atan(x / sqrt(1 - x * x));
	}


	/// Compute the arccosine
	inline real acos(real x) {

		if(abs(x) > 1) {
			UMATH_ERROR("acos", x, OUT_OF_DOMAIN);
			return nan();
		}

		if(x < 0)
			return atan(sqrt(1 - x * x) / x) + PI;
		else
			return atan(sqrt(1 - x * x) / x);
	}


	/// Compute the 2 argument arctangent
	inline real atan2(real y, real x) {

		if(x == 0) {

			if(y == 0) {
				UMATH_ERROR("atan2", y, OUT_OF_DOMAIN);
				return nan();
			}

			return sgn(y) * PI2;
		}

		if(x > 0) {
			return sgn(y) * atan(y / x);
		} else {
			return sgn(y) * atan(y / -x) + PI2;
		}
	}


	/// Compute the hyperbolic sine
	inline real sinh(real x) {
		real exp_x = exp(x);
		return (exp_x - 1.0 / exp_x) / 2.0;
	}


	/// Compute the hyperbolic cosine
	inline real cosh(real x) {
		real exp_x = exp(x);
		return (exp_x + 1.0 / exp_x) / 2.0;
	}


	// Compute the hyperbolic tangent
	inline real tanh(real x) {
		real exp_x = exp(x);
		return (exp_x - 1.0 / exp_x) / (exp_x + 1.0 / exp_x);
	}


	/// Compute the hyperbolic cotangent
	inline real coth(real x) {
		real exp_x = exp(x);
		return (exp_x + 1.0 / exp_x) / (exp_x - 1.0 / exp_x);
	}


	/// Compute the binomial coefficient
	inline long long binomial_coeff(unsigned int n, unsigned int m) {

		if(n < m) {
			UMATH_ERROR("binomial_coeff", n, IMPOSSIBLE_OPERATION);
			return 0;
		}

		// TO-DO Check out of range

		long long res = 1;

		for (int i = n; i > m; --i)
			res *= i;

		return res / fact(n - m);
	}


	/// Convert degrees to radians
	inline real radians(real degrees) {
		return degrees * DEG2RAD;
	}


	// Convert radians to degrees
	inline real degrees(real radians) {
		return radians * RAD2DEG;
	}


	/// Kronecker delta, equals 1 if i is equal to j, 0 otherwise
	template<typename T>
	inline T kronecker_delta(T i, T j, T tol = 0) {
		return (i - j < tol) ? 1 : 0;
	}

}

#endif
