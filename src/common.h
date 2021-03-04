#ifndef UROBORO_COMMON_H
#define UROBORO_COMMON_H

#include "./constants.h"
#include <cstdint>

namespace uroboro {

	// Calculate the square root of x using x86 Assembly instructions
	inline real sqrt(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			fsqrt
		}
		#else
		asm("fsqrt" : "+t"(x));
		return x;
		#endif
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
		return x;
		#endif
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
		return x;
		#endif
	}


	// Calculate tangent of x
	inline real tan(real x) {

		return sin(x) / cos(x);
	}


	// Calculate the cotangent of x
	inline real cot(real x) {

		return cos(x) / sin(x);
	}


	// Calculate the absolute value of x with x86 Assembly instructions
	inline real abs(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			fabs
		}
		#else
		asm("fabs" : "+t"(x));
		return x;
		#endif
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
