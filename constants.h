#ifndef UROBORO_CONSTANTS_H
#define UROBORO_CONSTANTS_H

#ifdef UROBORO_DOUBLE_PRECISION
using real = double;
#else
using real = float;
#endif

namespace uroboro {

	// PI constant
	constexpr real PI = 3.141593f;

	// PI constant
	constexpr double PI_PREC = 3.141592653589793238462643;

	// Half PI
	constexpr real PI2 = 1.570796f;

	// A quarter of PI.
	constexpr real PI4 = 0.785398f;

	// PI multiplied by 2.
	constexpr real PIDOUBLE = 6.283185f;

	// Inverse of the PI constant.
	constexpr real INVPI = 0.318310f;

	// Square root of PI constant.
	constexpr real SQRTPI = 1.772454f;

	// Euler constant.
	constexpr real E = 2.718282f;

	// Logarithms 2 of Euler constant.
	constexpr real LOG2E = 1.442695f;

	// Logarithm 10 of Euler constant.
	constexpr real LOG10E = 0.434295f;

	// Natural logarithm of 2.
	constexpr real LN2 = 0.693147f;

	// Natural logarithm of 10.
	constexpr real LN10 = 2.302585f;

	// Multiply a number by this constant to change it to radians.
	constexpr real DEG2RAD = 0.0174533f;

	// Multiply a number by this constant to change it to degrees.
	constexpr real RAD2DEG = 57.2958f;

	// Square root of 2.
	constexpr real SQRT2 = 1.414213f;

	// Inverse square root of 2.
	constexpr real INVSQR2 = 0.707106f;

}

#endif
