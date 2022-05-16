
///
/// @file core/constants.h Mathematical constants and algorithm parameters
///

#ifndef UROBORO_CONSTANTS_H
#define UROBORO_CONSTANTS_H


#include <limits>
#include <cstdint>


/// @macro UROBORO_X86 Define this macro if the
/// machine has a x86 architecture to use hardware
/// enhanced functions.


/// Order of Taylor series approximations
#ifndef UROBORO_TAYLOR_ORDER
#define UROBORO_TAYLOR_ORDER 9
#endif

/// Relative precision for derivative approximation
#ifndef UROBORO_DERIV_PREC
#define UROBORO_DERIV_PREC 100000.0
#endif

/// Default number of steps for integral approximation
#ifndef UROBORO_INTEGRATION_STEPS
#define UROBORO_INTEGRATION_STEPS 100
#endif

/// Biggest fractional part to ignore in powf computation
#ifndef UROBORO_POWF_APPROX_TOL
#define UROBORO_POWF_APPROX_TOL 0.00000001
#endif
	
/// Approximation tolerance for root finding
#ifndef UROBORO_ROOT_APPROX_TOL
#define UROBORO_ROOT_APPROX_TOL 0.00000001
#endif

/// Approximation tolerance for Newton's method root finding
#ifndef UROBORO_NEWTON_RAPHSON_TOL
#define UROBORO_NEWTON_RAPHSON_TOL 0.00000001
#endif

/// Approximation tolerance for bisection root finding
#ifndef UROBORO_BISECTION_APPROX_TOL
#define UROBORO_BISECTION_APPROX_TOL 0.00000001
#endif

/// Maximum number of iterations for bisection
#ifndef UROBORO_MAX_BISECTION_ITER
#define UROBORO_MAX_BISECTION_ITER 100
#endif

/// Maximum number of iterations for golden section search
#ifndef UROBORO_MAX_GOLDENSECTION_ITER
#define UROBORO_MAX_GOLDENSECTION_ITER 100
#endif

/// Maximum number of iterations for Halley's method
#ifndef UROBORO_MAX_HALLEY_ITER
#define UROBORO_MAX_HALLEY_ITER 100
#endif

/// Maximum number of iterations for Newton-Raphson
#ifndef UROBORO_MAX_NEWTON_ITER
#define UROBORO_MAX_NEWTON_ITER 100
#endif

/// Maximum number of iterations for Steffensen
#ifndef UROBORO_MAX_STEFFENSEN_ITER
#define UROBORO_MAX_STEFFENSEN_ITER 100
#endif

/// Maximum number of iterations for Chebyshev
#ifndef UROBORO_MAX_CHEBYSHEV_ITER
#define UROBORO_MAX_CHEBYSHEV_ITER 100
#endif

/// Maximum number of failed iterations for the Try-and-Catch algorithm
#ifndef UROBORO_MAX_TRYANDCATCH_ITER
#define UROBORO_MAX_TRYANDCATCH_ITER 100
#endif


#ifndef UROBORO_RAND_REAL_PREC
#define UROBORO_RAND_REAL_PREC (uint64_t(1) << 32)
#endif


/// @namespace uroboro Main namespace of the library which contains all functions and objects
namespace uroboro {


	/// A real number, defined as a floating point type.
	///
	/// The underlying type is determined by the defined macros:
	/// By default, `real` will be defined as the `double` type.
	/// If `UROBORO_FLOAT_PREC` is defined, `real` will be defined as a `float`,
	/// if `UROBORO_LONG_DOUBLE_PREC` is defined, `real` will be defined as a `long double`
	/// @note The `UROBORO_ARBITRARY_PREC` option is currently unsupported

#ifdef UROBORO_LONG_DOUBLE_PREC

	using real = long double;

#elif defined(UROBORO_FLOAT_PREC)

	using real = float;

#elif defined(UROBORO_ARBITRARY_PREC)

// TO-DO bigfloat arbitrary precision

#else

	using real = double;

#endif


	/// Machine epsilon for the real type
	constexpr real MACH_EPSILON = std::numeric_limits<real>::epsilon();

	/// The Phi (Golden Section) mathematical constant
	constexpr real PHI = 1.6180339887498948482045868;

	/// The inverse of the Golden Section mathematical constant
	constexpr real INVPHI = 0.6180339887498948482045868;

	/// The Pi mathematical constant
	constexpr real PI = 3.141592653589793238462643;

	/// Half of Pi
	constexpr real PI2 = PI / 2.0;

	/// A quarter of Pi
	constexpr real PI4 = PI / 4.0;

	/// Pi multiplied by 2
	constexpr real PIDOUBLE = PI * 2;

	/// The Tau mathematical constant (Pi times 2)
	constexpr real TAU = PI * 2;

	/// The inverse of Pi
	constexpr real INVPI = 1.0 / PI;

	/// The square root of Pi
	constexpr real SQRTPI = 1.772454;

	/// The Euler mathematical constant (e)
	constexpr real E = 2.718281828459045235360287;

	/// The binary logarithm of e
	constexpr real LOG2E = 1.44269504088896338700465094;

	/// The binary logarithm of 10
	constexpr real LOG210 = 3.32192809488736218170856773213;

	/// The base-10 logarithm of e
	constexpr real LOG10E = 0.434294481903;

	/// The natural logarithm of 2
	constexpr real LN2 = 0.69314718056;

	/// The natural logarithm of 10
	constexpr real LN10 = 2.30258509299;

	/// The scalar conversion factor from degrees to radians
	constexpr real DEG2RAD = 0.017453292519943295474371680598;

	/// The scalar conversion factor from radians to degrees
	constexpr real RAD2DEG = 57.2957795130823228646477218717;

	/// The square root of 2
	constexpr real SQRT2 = 1.41421356237;

	/// The inverse of the square root of 2
	constexpr real INVSQR2 = 0.707106781187;

	/// Order of Taylor series approximations
	constexpr int TAYLOR_ORDER = UROBORO_TAYLOR_ORDER;

	/// Default number of steps for integral approximation
	constexpr int INTEGRATION_STEPS = UROBORO_INTEGRATION_STEPS;

	/// Relative precision for derivative approximation
	constexpr real DERIV_PREC = UROBORO_DERIV_PREC;

	/// Biggest fractional part to ignore in powf computation
	constexpr real POWF_POWER_TOLERANCE = UROBORO_POWF_APPROX_TOL;

	/// Approximation tolerance for root finding
	constexpr real ROOT_APPROX_TOL = UROBORO_ROOT_APPROX_TOL;

	/// Approximation tolerance for the bisection algorithm
	constexpr real BISECTION_APPROX_TOL = UROBORO_BISECTION_APPROX_TOL;

	/// Approximation tolerance for the Newton-Raphson algorithm
	constexpr real NEWTON_RAPHSON_TOL = UROBORO_NEWTON_RAPHSON_TOL;

	/// Maximum number of iterations for the bisection algorithm
	constexpr unsigned int MAX_BISECTION_ITER = UROBORO_MAX_BISECTION_ITER;

	/// Maximum number of iterations for the golden section search algorithm
	constexpr unsigned int MAX_GOLDENSECTION_ITER = UROBORO_MAX_GOLDENSECTION_ITER;

	/// Maximum number of iterations for Halley's method
	constexpr unsigned int MAX_HALLEY_ITER = UROBORO_MAX_HALLEY_ITER;

	/// Maximum number of iterations for the Newton-Raphson algorithm
	constexpr unsigned int MAX_NEWTON_ITER = UROBORO_MAX_NEWTON_ITER;

	/// Maximum number of iterations for the Steffensen algorithm
	constexpr unsigned int MAX_STEFFENSEN_ITER = UROBORO_MAX_STEFFENSEN_ITER;

	/// Maximum number of iterations for the Chebyshev algorithm
	constexpr unsigned int MAX_CHEBYSHEV_ITER = UROBORO_MAX_CHEBYSHEV_ITER;

	/// Maximum number of failed iterations for the Try-and-Catch algorithm
	constexpr unsigned int MAX_TRYANDCATCH_ITER = UROBORO_MAX_TRYANDCATCH_ITER;


	/// Default precision for rand_real()
	constexpr uint64_t RAND_REAL_PREC = UROBORO_RAND_REAL_PREC;

}

#endif
