#ifndef UROBORO_COMMON_H
#define UROBORO_COMMON_H

#include "constants.h"
#include <cstdint>

namespace uroboro {


	// Calculate the square root of x using x86 Assembly instructions
	inline real sqrt(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			fsqrt
			fst x
		}
		#else
		asm("fsqrt" : "+t"(x));
		#endif

		return x;
	}


	// Approximate e^x using power series
	inline real exp(real x) {

		real x2 = x * x;
		real x4 = x2 * x2;

		return 1.f + x
		+ (x2 / 2.f)
		+ (x2 * x / 6.f)
		+ (x4 / 24.f)
		+ (x4 * x / 120.f)
		+ (x4 * x2 / 720.f)
		+ (x4 * x2 * x / 5040.f)
		+ (x4 * x4 / 40320.f);
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
			return x;

		return res;
	}


	// Calculate the natural logarithm of x
	// NOTE: Not so precise for numbers bigger than 20
	// Use log(real, unsigned int) for better precision
	inline real log(real x) {

		real precalc_x = (x - 1) / (x + 1);
		real x2 = precalc_x * precalc_x;
		real x3 = precalc_x * precalc_x * precalc_x;
		real x4 = x2 * x2;
		real x8 = x4 * x4;

		real result = precalc_x +
						(x3 / 3.f) +
						(x3 * x2 / 5.f) +
						(x3 * x4 / 7.f) +
						(x8 * precalc_x / 9.f) +
						(x8 * x2 * precalc_x / 11.f) +
						(x8 * x4 * precalc_x / 13.f) +
						(x8 * x4 * x3 / 15.f) +
						(x8 * x8 * precalc_x / 17.f) +
						(x8 * x8 * x3 / 19.f) +
						(x8 * x8 * x4 * precalc_x / 21.f);

		return result * 2;
	}

	// Calculate the binary logarithm of x
	inline real log2(real x) {

		return log(x) / log(2);
	}


	// Calculate the common logarithm of x
	inline real log10(real x) {

		return log(x) / log(10);
	}


	// Calculate the natural logarithm of x
	inline real log(real x, unsigned int iterations) {

		real precalc_x = (x - 1) / (x + 1);
		real result = 0;

		for (int i = 0; i < 7; ++i)
			result += pow(precalc_x, (2 * i) + 1) / (i * 2.f + 1);

		return result * 2;
	}


	// Calculate the binary logarithm of x
	inline real log2(real x, unsigned int iterations) {

		return log(x, iterations) / 0.69314718056f;
	}


	// Calculate the common logarithm of x
	inline real log10(real x, unsigned int iterations) {

		// This function doesn't use a constant for log(10)
		// because it wasn't precise enough
		return log(x, iterations) / log(10);
	}


	// Calculate x^n (where n has decimals)
	// WARNING: This function is considerably slower than pow(int, int)
	// WARNING: This function currently works only for inputs lower than 20
	inline real powf(real x, real n) {

		return exp(n * log(x));
	}


	// Calculate x! using iteration
	inline unsigned int fact(unsigned int x) {

		unsigned int res = 1;
		for (int i = x; i > 1; --i)
			res *= i;

		return res;
	}


	// Calculate sin(x) using x86 Assembly instructions
	inline real sin(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			fsin
			fst x
		}
		#else
		asm("fsin" : "+t"(x));
		#endif

		return x;
	}


	// Calculate cos(x) using
	inline real cos(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			fcos
			fst x
		}
		#else
		asm("fcos" : "+t"(x));
		#endif

		return x;
	}


	// Calculate tangent of x
	inline real tan(real x) {

		return sin(x) / cos(x);
	}


	// Calculate the cotangent of x
	inline real cot(real x) {

		return cos(x) / sin(x);
	}


	// Calculate asin(x) with x in range -1/+1 using Taylor series
	inline real asin(real x) {
#ifdef UROBORO_DOUBLE_PRECISION
		return -(-0.69813170079773212 * x * x - 0.87266462599716477) * x;
#else
		return -(-0.6981317f * x * x - 0.8726646f) * x;
#endif
	}


	// Calculate acos(x) with x in range -1/+1
	inline real acos(real x) {
#ifdef UROBORO_DOUBLE_PRECISION
		return (-0.69813170079773212 * x * x - 0.87266462599716477) * x + 1.5707963267948966;
#else
		return (-0.69813170079773212 * x * x - 0.87266462599716477) * x + PI2;
#endif
	}


	// Calculate atan(x) with x in range -1/+1 using Taylor series
	inline real atan(real x) {

		real x2 = x * x;
		real x4 = x2 * x2;

		return x
		- (x2 * x / 3.f)
		+ (x4 * x / 5.f)
		- (x4 * x2 * x / 7.f)
		+ (x4 * x4 * x2 * x / 11.f);
	}


	// Calculate the absolute value of x with x86 Assembly instructions
	inline real abs(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			fabs
			fst x
		}
		#else
		asm("fabs" : "+t"(x));
		#endif

		return x;
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


	// Return x + y * z
	// NOTE: This functions uses assembly instructions to speed up calculations
	// inline real fma(real x, real y, real z) {
	//
	// }


}

#endif
