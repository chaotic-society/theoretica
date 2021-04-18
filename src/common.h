#ifndef UROBORO_COMMON_H
#define UROBORO_COMMON_H

#include "./constants.h"
#include <cstdint>

namespace uroboro {


	//Calculate x^2
	inline real square(real x) {
		return x * x;
	}


	//Calculate x^3
	inline real cube(real x) {
		return x * x * x;
	}


	// Calculate the square root of x using x86 Assembly instructions
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


	// Calculate sin(x) using x86 Assembly instructions
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

	// Calculate cos(x) using x86 Assembly instructions
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


	// Calculate tangent of x
	inline real tan(real x) {

		// fptan usage TO-DO

		real s, c;

		#ifdef MSVC_ASM

		// TO-DO

		#else
		asm ("fsincos" : "=t"(c), "=u"(s) : "0"(x));
		#endif

		return s / c;
	}


	// Calculate the cotangent of x
	inline real cot(real x) {

		real s, c;

		#ifdef MSVC_ASM

		// TO-DO

		#else
		asm ("fsincos" : "=t"(c), "=u"(s) : "0"(x));
		#endif

		return c / s;
	}


	// Calculate the absolute value of x using x86 Assembly instructions
	inline real abs(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			fabs
		}
		#else
		asm("fabs" : "+t"(x));
		#endif

		return x;
	}


	// Calculate y * log2(x) using x86 Assembly instructions
	inline real fyl2x(real x, real y) {

		real res;

		#ifdef MSVC_ASM

		// TO-DO

		#else
		asm ("fyl2x" : "=t"(res) : "0"(x), "u"(y) : "st(1)");
		#endif

		return res;
	}


	// Calculate 2^x - 1 using x86 Assembly instructions
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


	// Calculate log2(x) using x86 Assembly instructions
	inline real log2(real x) {
		return fyl2x(x, 1.0);
	}


	// Calculate log10(x) using x86 Assembly instructions
	inline real log10(real x) {
		return fyl2x(x, 1.0 / LOG210);
	}


	// Calculate ln(x) using x86 Assembly instructions
	inline real ln(real x) {
		return fyl2x(x, 1.0 / LOG2E);
	}


	// Calculate x^n (where n is natural) using iteration
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


	// Clamp x between a and b
	inline real clamp(real x, real a, real b) {
		return x > b ? b : (x < a ? a : x);
	}


	// Return the biggest number between x and y
	inline real max(real x, real y) {
		return x > y ? x : y;
	}


	// Return the smallest number between x and y
	inline real min(real x, real y) {
		return x > y ? y : x;
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
