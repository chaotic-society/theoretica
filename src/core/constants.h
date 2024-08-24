
///
/// @file constants.h Mathematical constants and default algorithm
/// parameters. You may change the library's default behavior
/// by redefining the macros starting with THEORETICA_
///

#ifndef THEORETICA_CONSTANTS_H
#define THEORETICA_CONSTANTS_H


#include <limits>
#include <cstdint>

/// THEORETICA_DISABLE_X86 Define this macro to disable
/// Assembly x86 optimizations.
#ifndef THEORETICA_DISABLE_X86

#ifndef THEORETICA_X86
#if defined(__amd64__) || defined(__amd64) || defined(__x86_64__) \
	|| defined(_M_AMD64) || defined(_M_X64) || defined(i386) \
	|| defined(__i386) || defined(__i386__) || defined(_M_IX86) \
	|| defined(_X86_) || defined(_M_I86) || defined(__X86__)

/// THEORETICA_X86 This macro is automatically defined
/// on most compilers to enable Assembly x86 optimizations.
/// @see THEORETICA_DISABLE_X86
#define THEORETICA_X86
#endif
#endif
#endif


/// Tolerance for checking the elements of a
/// matrix (such as algebra::is_square)
#ifndef THEORETICA_ALGEBRA_ELEMENT_TOL
#define THEORETICA_ALGEBRA_ELEMENT_TOL (10*MACH_EPSILON)
#endif


/// Tolerance for eigensolvers
#ifndef THEORETICA_ALGEBRA_EIGEN_TOL
#define THEORETICA_ALGEBRA_EIGEN_TOL 1E-08
#endif


/// Maximum number of iterations for eigensolvers
#ifndef THEORETICA_ALGEBRA_EIGEN_ITER
#define THEORETICA_ALGEBRA_EIGEN_ITER 100
#endif


/// Order of Taylor series approximations
#ifndef THEORETICA_CORE_TAYLOR_ORDER
#define THEORETICA_CORE_TAYLOR_ORDER 12
#endif

/// Default number of steps for integral approximation
#ifndef THEORETICA_CALCULUS_INTEGRAL_STEPS
#define THEORETICA_CALCULUS_INTEGRAL_STEPS 100
#endif

// Default tolerance for integral approximation
#ifndef THEORETICA_CALCULUS_INTEGRAL_TOL
#define THEORETICA_CALCULUS_INTEGRAL_TOL 1E-08
#endif
	
/// Approximation tolerance for root finding
#ifndef THEORETICA_OPTIMIZATION_TOL
#define THEORETICA_OPTIMIZATION_TOL 1E-08
#endif

/// Maximum number of iterations for bisection
#ifndef THEORETICA_OPTIMIZATION_BISECTION_ITER
#define THEORETICA_OPTIMIZATION_BISECTION_ITER 100
#endif

/// Maximum number of iterations for golden section search
#ifndef THEORETICA_OPTIMIZATION_GOLDENSECTION_ITER
#define THEORETICA_OPTIMIZATION_GOLDENSECTION_ITER 100
#endif

/// Maximum number of iterations for Halley's method
#ifndef THEORETICA_OPTIMIZATION_HALLEY_ITER
#define THEORETICA_OPTIMIZATION_HALLEY_ITER 100
#endif

/// Maximum number of iterations for Newton-Raphson root finding
#ifndef THEORETICA_OPTIMIZATION_NEWTON_ITER
#define THEORETICA_OPTIMIZATION_NEWTON_ITER 100
#endif

/// Maximum number of iterations for Steffensen root finding
#ifndef THEORETICA_OPTIMIZATION_STEFFENSEN_ITER
#define THEORETICA_OPTIMIZATION_STEFFENSEN_ITER 100
#endif

/// Maximum number of iterations for Chebyshev root finding
#ifndef THEORETICA_OPTIMIZATION_CHEBYSHEV_ITER
#define THEORETICA_OPTIMIZATION_CHEBYSHEV_ITER 100
#endif

/// Maximum number of failed iterations for the Try-and-Catch algorithm
#ifndef THEORETICA_STATISTICS_TRYANDCATCH_ITER
#define THEORETICA_STATISTICS_TRYANDCATCH_ITER 100
#endif

/// Default variation for derivative approximation
#ifndef THEORETICA_CALCULUS_DERIV_STEP
#define THEORETICA_CALCULUS_DERIV_STEP 1E-3
#endif

/// Default step size for gradient descent minimization
#ifndef THEORETICA_OPTIMIZATION_MINGRAD_GAMMA
#define THEORETICA_OPTIMIZATION_MINGRAD_GAMMA -0.005
#endif

/// Default tolerance for gradient descent minimization
#ifndef THEORETICA_OPTIMIZATION_MINGRAD_TOLERANCE
#define THEORETICA_OPTIMIZATION_MINGRAD_TOLERANCE 1E-3
#endif

/// Maximum number of iterations for gradient descent minimization
#ifndef THEORETICA_OPTIMIZATION_MINGRAD_ITER
#define THEORETICA_OPTIMIZATION_MINGRAD_ITER 50000
#endif


#ifndef THEORETICA_STATISTICS_RAND_PREC

/// Default precision for random number generation using rand_uniform()
#ifdef THEORETICA_FLOAT_PREC
#define THEORETICA_STATISTICS_RAND_PREC (uint64_t(1) << 23)
#else
#define THEORETICA_STATISTICS_RAND_PREC (uint64_t(1) << 31)
#endif

#endif


/// Default depth of the Metropolis algorithm
#ifndef THEORETICA_STATISTICS_METROPOLIS_DEPTH
#define THEORETICA_STATISTICS_METROPOLIS_DEPTH 16
#endif


/// Enable constexpr in function declarations if C++14 is supported.
#if (__cplusplus >= 201402L)
#define TH_CONSTEXPR constexpr
#else
#define TH_CONSTEXPR
#endif


/// Enable constexpr in if statements if C++17 is supported.
#if (__cplusplus >= 201703L)
#define TH_CONSTIF constexpr
#else
#define TH_CONSTIF
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


	// Mathematical constants and default algorithm parameters.


	/// Machine epsilon for the real type
	constexpr real MACH_EPSILON = std::numeric_limits<real>::epsilon();

	/// The Phi (Golden Section) mathematical constant
	constexpr real PHI = 1.6180339887498948482045868;

	/// The inverse of the Golden Section mathematical constant
	constexpr real INVPHI = 0.6180339887498948482045868;

	/// The Pi mathematical constant
	constexpr real PI = 3.141592653589793238462643;

	/// Half of Pi
	constexpr real PI2 = 1.57079632679489655799898;

	/// A quarter of Pi
	constexpr real PI4 = PI / 4.0;

	/// Pi multiplied by 2
	constexpr real PIDOUBLE = PI * 2;

	/// The Tau mathematical constant (Pi times 2)
	constexpr real TAU = PI * 2;

	/// The inverse of Pi
	constexpr real INVPI = 1.0 / PI;

	/// The square root of Pi
	constexpr real SQRTPI = 1.7724538509055159927;

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
	constexpr real SQRT2 = 1.4142135623730950488;

	/// The inverse of the square root of 2
	constexpr real INVSQR2 = 0.7071067811865475;

	/// The square root of 3
	constexpr real SQRT3 = 1.732050807568877;

	/// Tolerance for the elements of matrices
	constexpr real ALGEBRA_ELEMENT_TOL = THEORETICA_ALGEBRA_ELEMENT_TOL;

	/// Tolerance for eigensolvers
	constexpr real ALGEBRA_EIGEN_TOL = THEORETICA_ALGEBRA_EIGEN_TOL;

	/// Maximum number of iterations for eigensolvers
	constexpr real ALGEBRA_EIGEN_ITER = THEORETICA_ALGEBRA_EIGEN_ITER;

	/// Order of Taylor series approximations
	constexpr int CORE_TAYLOR_ORDER = THEORETICA_CORE_TAYLOR_ORDER;

	/// Default number of steps for integral approximation
	constexpr int CALCULUS_INTEGRAL_STEPS = THEORETICA_CALCULUS_INTEGRAL_STEPS;

	// Default tolerance for integral approximation
	constexpr real CALCULUS_INTEGRAL_TOL = THEORETICA_CALCULUS_INTEGRAL_TOL;

	/// Approximation tolerance for root finding
	constexpr real OPTIMIZATION_TOL = THEORETICA_OPTIMIZATION_TOL;

	/// Maximum number of iterations for the bisection algorithm
	constexpr unsigned int OPTIMIZATION_BISECTION_ITER = THEORETICA_OPTIMIZATION_BISECTION_ITER;

	/// Maximum number of iterations for the golden section search algorithm
	constexpr unsigned int OPTIMIZATION_GOLDENSECTION_ITER = THEORETICA_OPTIMIZATION_GOLDENSECTION_ITER;

	/// Maximum number of iterations for Halley's method
	constexpr unsigned int OPTIMIZATION_HALLEY_ITER = THEORETICA_OPTIMIZATION_HALLEY_ITER;

	/// Maximum number of iterations for the Newton-Raphson algorithm
	constexpr unsigned int OPTIMIZATION_NEWTON_ITER = THEORETICA_OPTIMIZATION_NEWTON_ITER;

	/// Maximum number of iterations for the Steffensen algorithm
	constexpr unsigned int OPTIMIZATION_STEFFENSEN_ITER = THEORETICA_OPTIMIZATION_STEFFENSEN_ITER;

	/// Maximum number of iterations for the Chebyshev algorithm
	constexpr unsigned int OPTIMIZATION_CHEBYSHEV_ITER = THEORETICA_OPTIMIZATION_CHEBYSHEV_ITER;

	/// Maximum number of failed iterations for the Try-and-Catch algorithm
	constexpr unsigned int STATISTICS_TRYANDCATCH_ITER = THEORETICA_STATISTICS_TRYANDCATCH_ITER;

	/// Default variation for derivative approximation
	constexpr real CALCULUS_DERIV_STEP = THEORETICA_CALCULUS_DERIV_STEP;

	/// Default step size for gradient descent minimization
	constexpr real OPTIMIZATION_MINGRAD_GAMMA = THEORETICA_OPTIMIZATION_MINGRAD_GAMMA;

	/// Default tolerance for gradient descent minimization
	constexpr real OPTIMIZATION_MINGRAD_TOLERANCE = THEORETICA_OPTIMIZATION_MINGRAD_TOLERANCE;

	/// Maximum number of iterations for gradient descent minimization
	constexpr unsigned int OPTIMIZATION_MINGRAD_ITER = THEORETICA_OPTIMIZATION_MINGRAD_ITER;

	/// Default precision for random number generation using rand_uniform()
	constexpr uint64_t STATISTICS_RAND_PREC = THEORETICA_STATISTICS_RAND_PREC;

	/// Default depth of the Metropolis algorithm
	constexpr unsigned int STATISTICS_METROPOLIS_DEPTH = THEORETICA_STATISTICS_METROPOLIS_DEPTH;

}

// Define THEORETICA_NO_NAMESPACE_ALIAS to prevent
// defining the alias "th" for "theoretica"
#ifndef THEORETICA_NO_NAMESPACE_ALIAS
/// @namespace th Alias for the theoretica namespace
namespace th = theoretica;
#endif

#endif
