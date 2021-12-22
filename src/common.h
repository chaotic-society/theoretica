#ifndef UROBORO_COMMON_H
#define UROBORO_COMMON_H

#include "./constants.h"


namespace uroboro {


	//Compute x^2
	inline real square(real x) {
		return x * x;
	}


	//Compute x^3
	inline real cube(real x) {
		return x * x * x;
	}


	// Compute the square root of x using x86 Assembly instructions
	inline real sqrt(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			fsqrt
		}
		#else
		asm("fsqrt" : "+t"(x));
		#endif

		return x;
	}


	// Compute the absolute value of x using x86 Assembly instructions
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


	// Return the sign of x (1 if positive, -1 if negative, 0 if null)
	inline int sgn(real x) {
		return (x > 0) ? 1 : (x < 0 ? -1 : 0);
	}


	// Return the greatest number between x and y
	inline real max(real x, real y) {
		
		#ifdef UROBORO_BRANCHLESS
			return (x + y + abs(x - y)) / 2.0;
		#else
			return x > y ? x : y;
		#endif
	}


	// Return the smallest number between x and y
	inline real min(real x, real y) {

		#ifdef UROBORO_BRANCHLESS
			return (x + y - abs(x - y)) / 2.0;
		#else
			return x > y ? y : x;
		#endif
	}


	// Clamp x between a and b
	inline real clamp(real x, real a, real b) {

#ifdef UROBORO_FORCE_BRANCHLESS

		// The branchless version might be slower or equal
		// on most compilers
		return min(max(x, a), b);
#else
		return x > b ? b : (x < a ? a : x);
#endif
	}


	// x86 instruction wrappers

#ifdef UROBORO_X86

	// Compute y * log2(x) using x86 Assembly instructions
	inline real fyl2x(real x, real y) {

		real res;

		#ifdef MSVC_ASM

		// TO-DO

		#else
		asm ("fyl2x" : "=t"(res) : "0"(x), "u"(y) : "st(1)");
		#endif

		return res;
	}


	// Compute 2^x - 1 using x86 Assembly instructions
	// Works only between -1 and +1
	// May become particularly incorrect near boundaries
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

	// Compute log2(x) using x86 Assembly instructions
	inline real log2(real x) {
		return fyl2x(x, 1.0);
	}


	// Compute log10(x) using x86 Assembly instructions
	inline real log10(real x) {
		return fyl2x(x, 1.0 / LOG210);
	}


	// Compute ln(x) using x86 Assembly instructions
	inline real ln(real x) {
		return fyl2x(x, 1.0 / LOG2E);
	}


	// Compute x^n (where n is natural) using iteration
	inline real pow(real x, int n) {

		real res;
		if(n > 0) {

			res = x;
			for(int i = 1; i < n; ++i)
				res *= x;

		} else if(n < 0) {

			res = 1 / x;
			for(int i = 1; i < -n; ++i)
				res /= x;
		} else
			return 1.0;

		return res;
	}


	constexpr real POW_APPROXIMATION_TOLERANCE = 0.000001;


	// Approximate e^x using x86 Assembly instructions
	// Works only for positive <x>
	inline real exp_approx(real x) {

		int x_int = uroboro::abs(int(x - 0.5));
		real x_fract = uroboro::abs(x - x_int);

		// Compute e^x as e^int(x) * e^fract(x)
		// Where e^fract(x) is calculated as 2^(fract(x) / ln2)
		return uroboro::pow(E, x_int) * square(f2xm1(x_fract / (2 * LN2)) + 1);
	}


	// Approximate powf(real, real) = x^a
	// Using pow(x, int(a)) * exp(fract(a) * ln(x))
	inline real powf_approx(real x, real a) {

		if(a < 0)
			return 1.0 / powf_approx(x, abs(a));

		int a_int = uroboro::abs(int(a - 0.5));
		real a_fract = uroboro::abs(a - a_int);
		real x_int_pwr = uroboro::pow(x, a_int);

		// Compute x^fract(a) as e^(x * log2(fract(a) / log2(e)))
		return x_int_pwr * (a_fract >= POW_APPROXIMATION_TOLERANCE ?
			exp_approx(fyl2x(x, a_fract / LOG2E)) : 1);
	}


	// Compute e^x
	inline real exp(real x) {
		return powf_approx(E, x);
	}


	// Compute sin(x) using x86 Assembly instructions
	inline real sin(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			fsin
		}
		#else
		asm("fsin" : "+t"(x));
		#endif

		return x;
	}

	// Compute cos(x) using x86 Assembly instructions
	inline real cos(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			fcos
		}
		#else
		asm("fcos" : "+t"(x));
		#endif

		return x;
	}


	// Compute the tangent of x
	inline real tan(real x) {

		real s, c;

		#ifdef MSVC_ASM

		// TO-DO

		#else
		asm ("fsincos" : "=t"(c), "=u"(s) : "0"(x));
		#endif

		return s / c;
	}


	// Compute the cotangent of x
	inline real cot(real x) {

		real s, c;

		#ifdef MSVC_ASM

		// TO-DO

		#else
		asm ("fsincos" : "=t"(c), "=u"(s) : "0"(x));
		#endif

		return c / s;
	}


	// Compute the arctangent
	inline real atan(real x) {

		if(abs(x) > 1)
			return atan(PI2 - atan(1 / x));

		if(abs(x) < 0.6) {

			// Use Taylor in [-0.6, 0.6]

			real x2 = x * x;
			real x4 = x2 * x2;

			return x
				- (x2 * x / 3.f)
				+ (x4 * x / 5.f)
				- (x4 * x2 * x / 7.f)
				+ (x4 * x4 * x / 9.f);

		} else {

			// Fast approximation of atan in [-1, 1]
			// More accurate than Taylor in [-1, -0.6] and [0.6, 1]
			return PI4 * x - x * (abs(x) - 1) * (0.2447 + 0.0663 * abs(x));
		}
	}


	// Compute the arcsine
	inline real asin(real x) {
		return atan(x / sqrt(1 - x * x));
	}


	// Compute the arccosine
	inline real acos(real x) {
		return atan(sqrt(1 - x * x) / x);
	}


	// Compute the 2 argument arctangent
	inline real atan2(real y, real x) {

		if(x == 0)
			return sgn(y) * PI2;
			// throw y == 0 ...

		return atan(y / x) - clamp(sgn(x), -1, 0) * PI * sgn(y);
	}


	// Compute the hyperbolic sine
	inline real sinh(real x) {
		return (exp(x) - exp(-x)) / 2.0;
	}


	// Compute the hyperbolic cosine
	inline real cosh(real x) {
		return (exp(x) + exp(-x)) / 2.0;
	}


	// Compute the hyperbolic tangent
	inline real tanh(real x) {
		real exp_2x = exp(2 * x);
		return (exp_2x - 1) / (exp_2x + 1);
	}


	// Compute the factorial of n
	inline long long fact(unsigned int n) {

		long long res = 1;
		for (int i = n; i > 1; --i)
			res *= i;

		return res;
	}


	// Compute the binomial coefficient
	inline long long binomial_coeff(unsigned int n, unsigned int m) {

		if(n < m)
			// throw...
			return 0;

		long long res = 1;

		for (int i = n; i > m; --i)
			res *= i;

		return res / fact(n - m);
	}


	// Convert degrees to radians
	inline real radians(real degrees) {
		return degrees * DEG2RAD;
	}


	// Convert radians to degrees
	inline real degrees(real radians) {
		return radians * RAD2DEG;
	}

}

#endif
