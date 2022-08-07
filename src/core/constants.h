
///
/// @file constants.h Mathematical constants and algorithm parameters
///

#ifndef THEORETICA_CONSTANTS_H
#define THEORETICA_CONSTANTS_H


#include <limits>
#include <cstdint>


/// @macro THEORETICA_X86 Define this macro if the
/// machine has a x86 architecture to use hardware
/// enhanced functions.
#ifndef THEORETICA_DISABLE_X86
#ifndef THEORETICA_X86
#if defined(__amd64__) || defined(__amd64) || defined(__x86_64__) \
	|| defined(_M_AMD64) || defined(_M_X64) || defined(i386) \
	|| defined(__i386) || defined(__i386__) || defined(_M_IX86) \
	|| defined(_X86_) || defined(_M_I86) || defined(__X86__)	
#define THEORETICA_X86
#endif
#endif
#endif


/// Order of Taylor series approximations
#ifndef THEORETICA_TAYLOR_ORDER
#define THEORETICA_TAYLOR_ORDER 12
#endif

/// Relative precision for derivative approximation
#ifndef THEORETICA_DERIV_PREC
#define THEORETICA_DERIV_PREC 100000.0
#endif

/// Default number of steps for integral approximation
#ifndef THEORETICA_INTEGRATION_STEPS
#define THEORETICA_INTEGRATION_STEPS 100
#endif

/// Biggest fractional part to ignore in powf computation
#ifndef THEORETICA_POWF_APPROX_TOL
#define THEORETICA_POWF_APPROX_TOL 0.00000001
#endif
	
/// Approximation tolerance for root finding
#ifndef THEORETICA_ROOT_APPROX_TOL
#define THEORETICA_ROOT_APPROX_TOL 0.00000001
#endif

/// Approximation tolerance for Newton's method root finding
#ifndef THEORETICA_NEWTON_RAPHSON_TOL
#define THEORETICA_NEWTON_RAPHSON_TOL 0.00000001
#endif

/// Approximation tolerance for bisection root finding
#ifndef THEORETICA_BISECTION_APPROX_TOL
#define THEORETICA_BISECTION_APPROX_TOL 0.00000001
#endif

/// Maximum number of iterations for bisection
#ifndef THEORETICA_MAX_BISECTION_ITER
#define THEORETICA_MAX_BISECTION_ITER 100
#endif

/// Maximum number of iterations for golden section search
#ifndef THEORETICA_MAX_GOLDENSECTION_ITER
#define THEORETICA_MAX_GOLDENSECTION_ITER 100
#endif

/// Maximum number of iterations for Halley's method
#ifndef THEORETICA_MAX_HALLEY_ITER
#define THEORETICA_MAX_HALLEY_ITER 100
#endif

/// Maximum number of iterations for Newton-Raphson
#ifndef THEORETICA_MAX_NEWTON_ITER
#define THEORETICA_MAX_NEWTON_ITER 100
#endif

/// Maximum number of iterations for Steffensen
#ifndef THEORETICA_MAX_STEFFENSEN_ITER
#define THEORETICA_MAX_STEFFENSEN_ITER 100
#endif

/// Maximum number of iterations for Chebyshev
#ifndef THEORETICA_MAX_CHEBYSHEV_ITER
#define THEORETICA_MAX_CHEBYSHEV_ITER 100
#endif

/// Maximum number of failed iterations for the Try-and-Catch algorithm
#ifndef THEORETICA_MAX_TRYANDCATCH_ITER
#define THEORETICA_MAX_TRYANDCATCH_ITER 100
#endif

/// Default variation for derivative approximation
#ifndef THEORETICA_DERIV_STEPSIZE
#define THEORETICA_DERIV_STEPSIZE 0.001
#endif

/// Default step size for gradient descent minimization
#ifndef THEORETICA_MINGRAD_GAMMA
#define THEORETICA_MINGRAD_GAMMA -0.005
#endif

/// Default tolerance for gradient descent minimization
#ifndef THEORETICA_MINGRAD_TOLERANCE
#define THEORETICA_MINGRAD_TOLERANCE 0.001
#endif

/// Maximum number of iterations for gradient descent minimization
#ifndef THEORETICA_MINGRAD_MAX_ITER
#define THEORETICA_MINGRAD_MAX_ITER 50000
#endif


#ifndef THEORETICA_RAND_REAL_PREC

/// Default precision for random number generation using rand_uniform()
#ifdef THEORETICA_FLOAT_PREC
#define THEORETICA_RAND_REAL_PREC (uint64_t(1) << 23)
#else
#define THEORETICA_RAND_REAL_PREC (uint64_t(1) << 32)
#endif

#endif


/// Enable constexpr in function declarations if C++14 is supported
#if (__cplusplus >= 201402L)
#define TH_CONSTEXPR constexpr
#else
#define TH_CONSTEXPR
#endif


/// @namespace theoretica Main namespace of the library which contains all functions and objects
namespace theoretica {


	/// A real number, defined as a floating point type.
	///
	/// The underlying type is determined by the defined macros:
	/// By default, `real` will be defined as the `double` type.
	/// If `THEORETICA_FLOAT_PREC` is defined, `real` will be defined as a `float`,
	/// if `THEORETICA_LONG_DOUBLE_PREC` is defined, `real` will be defined as a `long double`
	/// @note The `THEORETICA_ARBITRARY_PREC` option is currently unsupported

#ifdef THEORETICA_LONG_DOUBLE_PREC

	using real = long double;

#elif defined(THEORETICA_FLOAT_PREC)

	using real = float;

#elif defined(THEORETICA_ARBITRARY_PREC)

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
	constexpr real SQRT2 = 1.414213562373095;

	/// The inverse of the square root of 2
	constexpr real INVSQR2 = 0.7071067811865475;

	/// The square root of 3
	constexpr real SQRT3 = 1.732050807568877;

	/// Order of Taylor series approximations
	constexpr int TAYLOR_ORDER = THEORETICA_TAYLOR_ORDER;

	/// Default number of steps for integral approximation
	constexpr int INTEGRATION_STEPS = THEORETICA_INTEGRATION_STEPS;

	/// Relative precision for derivative approximation
	constexpr real DERIV_PREC = THEORETICA_DERIV_PREC;

	/// Biggest fractional part to ignore in powf computation
	constexpr real POWF_POWER_TOLERANCE = THEORETICA_POWF_APPROX_TOL;

	/// Approximation tolerance for root finding
	constexpr real ROOT_APPROX_TOL = THEORETICA_ROOT_APPROX_TOL;

	/// Approximation tolerance for the bisection algorithm
	constexpr real BISECTION_APPROX_TOL = THEORETICA_BISECTION_APPROX_TOL;

	/// Approximation tolerance for the Newton-Raphson algorithm
	constexpr real NEWTON_RAPHSON_TOL = THEORETICA_NEWTON_RAPHSON_TOL;

	/// Maximum number of iterations for the bisection algorithm
	constexpr unsigned int MAX_BISECTION_ITER = THEORETICA_MAX_BISECTION_ITER;

	/// Maximum number of iterations for the golden section search algorithm
	constexpr unsigned int MAX_GOLDENSECTION_ITER = THEORETICA_MAX_GOLDENSECTION_ITER;

	/// Maximum number of iterations for Halley's method
	constexpr unsigned int MAX_HALLEY_ITER = THEORETICA_MAX_HALLEY_ITER;

	/// Maximum number of iterations for the Newton-Raphson algorithm
	constexpr unsigned int MAX_NEWTON_ITER = THEORETICA_MAX_NEWTON_ITER;

	/// Maximum number of iterations for the Steffensen algorithm
	constexpr unsigned int MAX_STEFFENSEN_ITER = THEORETICA_MAX_STEFFENSEN_ITER;

	/// Maximum number of iterations for the Chebyshev algorithm
	constexpr unsigned int MAX_CHEBYSHEV_ITER = THEORETICA_MAX_CHEBYSHEV_ITER;

	/// Maximum number of failed iterations for the Try-and-Catch algorithm
	constexpr unsigned int MAX_TRYANDCATCH_ITER = THEORETICA_MAX_TRYANDCATCH_ITER;

	/// Default variation for derivative approximation
	constexpr real DERIV_STEPSIZE = THEORETICA_DERIV_STEPSIZE;

	/// Default step size for gradient descent minimization
	constexpr real MINGRAD_GAMMA = THEORETICA_MINGRAD_GAMMA;

	/// Default tolerance for gradient descent minimization
	constexpr real MINGRAD_TOLERANCE = THEORETICA_MINGRAD_TOLERANCE;

	/// Maximum number of iterations for gradient descent minimization
	constexpr unsigned int MINGRAD_MAX_ITER = THEORETICA_MINGRAD_MAX_ITER;

	/// Default precision for random number generation using rand_uniform()
	constexpr uint64_t RAND_REAL_PREC = THEORETICA_RAND_REAL_PREC;

}

#endif
